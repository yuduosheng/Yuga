#include "FULLIntegrator.h"
#include "TIMING_BREAKDOWN.h"
ofstream debug("debug.txt");

FullIntegrator::FullIntegrator(TetMesh *tetMesh, double dampingMassCoef, double dampingStiffnessCoef, int maxIterations, double timestep)
{
	_tetMesh = tetMesh;
	_material = new MaterialLinearElasticity(tetMesh);
	_massesMatrix = _tetMesh->getMassesMatrix();

	//vector<TRIPLET>& matrix = _massesMatrix.matrix();
	//cout << matrix.size() << " " << _massesMatrix.rows() * _massesMatrix.rows() << endl;
	//for (int i = 0; i < matrix.size(); i++)
	//{
	//	debug << matrix[i].row() << "  " << matrix[i].col() << "  " << matrix[i].value() << endl;
	//}
	//cout << endl;


	_vDofs = _tetMesh->GetVetexNumber() * 3;
	_dampingMassCoef = dampingMassCoef;
	_dampingStiffnessCoef = dampingStiffnessCoef;
	_maxIterations = maxIterations;
	_timestep = timestep;

	_externalForces.resize(_vDofs);
	_externalForces.setZero();
	_internalForces.resize(_vDofs);
	_internalForces.setZero();
	//cout << _internalForces.size() << endl;
	//COO_MATRIX stiffnessTest;
	//stiffnessTest.resize(12, 12);
	//_material->computeElementStiffnessMatrix(0, stiffnessTest);
	//
	//MATRIX& pFpu = _tetMesh->getTet(0)._PFPu;
	//debug << pFpu << endl;
	//debug << "=================" <<endl;
	//vector<TRIPLET>& matrix = stiffnessTest.matrix();
	//cout << matrix.size() << " " << _vDofs * _vDofs << endl;
	//for (int i = 0; i < matrix.size(); i++)
	//{
	//	debug << matrix[i].row() << "  " << matrix[i].col() << "  " << matrix[i].value() << endl;
	//}
	//cout << endl;


	_material->computeStiffnessMatrix(_stiffnessMatrix);
	//vector<TRIPLET>& matrix = _stiffnessMatrix.matrix();
	//cout << matrix.size() << " " << _vDofs * _vDofs << endl;
	//for (int i = 0; i < matrix.size(); i++)
	//{
	//	debug << matrix[i].row() << "  " << matrix[i].col() << "  " << matrix[i].value() << endl;
	//}
	//cout << endl;

	_epsilon = 0.1;

	q.resize(_vDofs);
	qvel.resize(_vDofs);
	q_old.resize(_vDofs);
	qvel_old.resize(_vDofs);
	qresidual.resize(_vDofs);
	deltaV.resize(_vDofs);

	_dampingMatrix = _massesMatrix * dampingMassCoef +  _stiffnessMatrix * dampingStiffnessCoef;
	//vector<TRIPLET>& matrix = _dampingMatrix.matrix();
	//cout << matrix.size() << " " << _dampingMatrix.rows() * _dampingMatrix.rows() << endl;
	//for (int i = 0; i < matrix.size(); i++)
	//{
	//	debug << matrix[i].row() << "  " << matrix[i].col() << "  " << matrix[i].value() << endl;
	//}
	//cout << endl;


	q.setZero();
	qvel.setZero();

	jacobiPreconditionedCGSolver = new CGSolver(true);
}

void FullIntegrator::SetupSystemMatrix(int numIter)
{
	_systemMatrix = _stiffnessMatrix;
	if (numIter != 0)
	{
		deltaV = q_old - q;

		_systemMatrix.multiplyVector(deltaV, qresidual);//qresidual = _stiffnessMatrix * (q_o - q)
	}
	// build effective stiffness: 
	// Keff = M + h D + h^2 * K
	// compute force residual, store it into aux variable qresidual
	// qresidual = h * (-D qdot - fint + fext - h * K * qdot)) // this is semi-implicit Euler
	// qresidual = M (qvel_1 - qvel) + h * (-D qdot - fint + fext - K * (q_1 - q + h qdot) )) // for fully implicit Euler

	_systemMatrix *= _timestep;
	_systemMatrix.plus(_dampingMatrix);
	_systemMatrix.multiplyVectorAdd(qvel, qresidual);//qresidual += (h*K + D)qdot
	_systemMatrix *= _timestep; // h^2 * K + h * D
	_systemMatrix.plus(_massesMatrix);

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
	//debug << _internalForces << endl;

	do
	{
		TIMING_BREAKDOWN::tic();
		SetupSystemMatrix(numIter);
		TIMING_BREAKDOWN::toc("Setup System Matrix");
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
	  TIMING_BREAKDOWN::tic();
	  //solve
	  jacobiPreconditionedCGSolver->solve(_systemMatrix ,qresidual, deltaV);
	  //update state
	  for(int i=0; i < _vDofs; i++)
      {
          qvel[i] += deltaV[i];
          q[i] = q_old[i] + _timestep * qvel[i];
      }	  
	  numIter++;
	  TIMING_BREAKDOWN::toc("Solve System");

	} while (numIter < _maxIterations);

	TIMING_BREAKDOWN::tic();
	_tetMesh->setCurPosition(q);
	TIMING_BREAKDOWN::toc("Set Position");
}

void FullIntegrator::addVertexForce(int index, VEC3F& force)
{
	_externalForces[index * 3] += force[0];
	_externalForces[index * 3 + 1] += force[1];
	_externalForces[index * 3 + 2] += force[2];
}