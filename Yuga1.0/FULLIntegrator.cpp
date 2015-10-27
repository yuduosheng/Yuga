#include "FULLIntegrator.h"

FullIntegrator::FullIntegrator(TetMesh *tetMesh, double dampingMassCoef, double dampingStiffnessCoef, int maxIterations, double timestep)
{
	_tetMesh = tetMesh;
	_material = new MaterialLinearElasticity(tetMesh, 1e7, 0.4);
	_massesMatrix = _tetMesh->getMassesMatrix();

	_vDofs = _tetMesh->GetVetexNumber * 3;
	_dampingMassCoef = dampingMassCoef;
	_dampingStiffnessCoef = dampingStiffnessCoef;
	_maxIterations = maxIterations;
	_timestep = timestep;

	_externalForces.resize(_vDofs);
	_externalForces.setZero();
	_internalForces.resize(_vDofs);
	_internalForces.setZero();

	_material->computeStiffnessMatrix(_stiffnessMatrix);
	_epsilon = 0.1;

	q.resize(_tetMesh->GetVetexNumber() * 3);
	qvel.resize(_tetMesh->GetVetexNumber() * 3);
	q_old.resize(_tetMesh->GetVetexNumber() * 3);
	qvel_old.resize(_tetMesh->GetVetexNumber() * 3);

	jacobiPreconditionedCGSolver = new CGSolver(_stiffnessMatrix);
}

void FullIntegrator::SetupSystemMatrix(int numIter)
{

	if (numIter != 0)
	{
		deltaV = q_old - q;

		_stiffnessMatrix.multiplyVector(deltaV, qresidual);//qresidual = _stiffnessMatrix * (q_o - q)
	}
	// build effective stiffness: 
	// Keff = M + h D + h^2 * K
	// compute force residual, store it into aux variable qresidual
	// qresidual = h * (-D qdot - fint + fext - h * K * qdot)) // this is semi-implicit Euler
	// qresidual = M (qvel_1 - qvel) + h * (-D qdot - fint + fext - K * (q_1 - q + h qdot) )) // for fully implicit Euler

	_stiffnessMatrix *= _timestep;
	_stiffnessMatrix.plus(_dampingMatrix);
	_stiffnessMatrix.multiplyVectorAdd(qvel, qresidual);//qresidual += (h*K + D)qdot
	_stiffnessMatrix *= _timestep; // h^2 * K + h * D
	_stiffnessMatrix.plus(_massesMatrix);

	// add externalForces, internalForces
	for (int i = 0; i < _vDofs; i++)
	{
		qresidual[i] += _internalForces[i] - _externalForces[i];
		qresidual[i] *= -_timestep;
	}

	if (numIter != 0) // can skip on first iteration (zero contribution)
	{
		// add M * (qvel_1 - qvel) to qresidual
		deltaV = qvel_old - qvel;

		_massesMatrix.multiplyVectorAdd(deltaV, qresidual);
	}
}

void FullIntegrator::DoTimeStep()
{
	int numIter = 0;
    double error0 = 0; // error after the first step
    double errorQuotient;

    // store current amplitudes and set initial guesses for qaccel, qvel
    q_old = q; 
    qvel_old = qvel;
	
	_material->computeInternalForce(_internalForces);

	do
	{
		SetupSystemMatrix(numIter);
		
		double error = 0.0;
   		for(int i = 0; i < _vDofs; i++)
      		error += qresidual[i] * qresidual[i];
			  
		    // on the first iteration, compute initial error
      if (numIter == 0) 
  	  {
     	 error0 = error;
     	 errorQuotient = 1.0;
      }
      else
   	  {
    	  // error divided by the initial error, before performing this iteration
    	  errorQuotient = error / error0; 
      }

      if (errorQuotient < _epsilon * _epsilon)
    	  break;
	  //solve
	  jacobiPreconditionedCGSolver->solvePCG(qresidual, deltaV);
	  //update state
	  for(int i=0; i < _vDofs; i++)
      {
          qvel[i] += deltaV[i];
          q[i] = q_old[i] + _timestep * qvel[i];
      }	  
	  numIter++;
	} while (numIter < _maxIterations);

	_tetMesh->setCurPosition(q);

}

void FullIntegrator::addVertexForce(int index, VEC3F& force)
{
	_externalForces[index * 3] += force[0];
	_externalForces[index * 3 + 1] += force[1];
	_externalForces[index * 3 + 2] += force[2];
}