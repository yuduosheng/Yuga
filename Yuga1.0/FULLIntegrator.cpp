#include "FULLIntegrator.h"

FullIntegrator::FullIntegrator(TetMesh *tetMesh, double dampingMassCoef = 0.25, double dampingStiffnessCoef = 0.5, int maxIterations = 5, double timestep)
{
	_tetMesh = tetMesh;
	_material = new MaterialLinearElasticity(1e7, 0.4);
	_internalForces = _material->getInternalForce();
	_stiffnessMatrix = _material->getStiffnessMatrix();
	_massesMatrix = _tetMesh->getMassesMatrix();

	_vDofs = _tetMesh->GetVetexNumber * 3;
	_dampingMassCoef = dampingMassCoef;
	_dampingStiffnessCoef = dampingStiffnessCoef;
	_maxIterations = maxIterations;
	_timestep = timestep;
	
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
	
	do
	{
		SetupSystemMatrix(numIter)
		
		double error = 0.0;
   		for(int i = 0; i < _vDofs; i++)
      		error += qresidual * qresidual[i];
			  
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

      if (errorQuotient < epsilon * epsilon)
    	  break;
	  //solve
	  jacobiPreconditionedCGSolver->solvePCG(qresidual, deltaV);
	  //update state
	  for(int i=0; i<r; i++)
      {
          qvel[i] += deltaV[i];
          q[i] = q_1[i] + _timestep * qvel[i];
      }	  
	  numIter++;
	}
   	while(numIter < _maxIterations)
}