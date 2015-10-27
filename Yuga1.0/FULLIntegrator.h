#ifndef FULLINTEGRATOR_H
#define FULLINTEGRATOR_H
#include "TetMesh.h"
#include "Material.h"
#include "COO_MATRIX.h"
#include "MaterialLinearElasticity.h"
#include "CONJUGATE_GRADIENT.h"
typedef CONJUGATE_GRADIENT CGSlover;
class FullIntegrator
{
protected:
	VECTOR q;//diceplacement
	VECTOR qvel;//current velocities of deformation amplitudes
	//VECTOR qaccel;//current acceleration
	VECTOR qresidual, deltaV;//aux integration variables
	VECTOR q_old;//variables in previous time-step
	VECTOR qvel_old;
	//VECTOR qaccel_old;


	VECTOR _internalForces;
	VECTOR _externalForces;

	COO_MATRIX _massesMatrix;
	COO_MATRIX _dampingMatrix;
	COO_MATRIX _stiffnessMatrix;
	
	CGSlover *jacobiPreconditionedCGSolver;
	//these two store the damping parameters
	double _dampingMassCoef;
	double _dampingStiffnessCoef;

	double _timestep;
	double _vDofs;

	int _maxIterations;
	double _epsilon;

	TetMesh *_tetMesh;
	Material *_material;
public:
	FullIntegrator(TetMesh *tetMesh, double dampingMassCoef, double dampingStiffnessCoef, int maxIterations, double timestep);
	~FullIntegrator(){}

	virtual void SetupSystemMatrix(int numIter);
	virtual void DoTimeStep();
};
#endif