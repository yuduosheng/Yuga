#ifndef FULLINTEGRATOR_H
#define FULLINTEGRATOR_H
#include "TetMesh.h"
#include "Material.h"
#include "COO_MATRIX.h"
#include "MaterialLinearElasticity.h"
class FullIntegrator
{
protected:
	VECTOR q;//diceplacement
	VECTOR qvel;//current velocities of deformation amplitudes
	VECTOR qaccel;//current acceleration
	VECTOR qresidul, qdelta;//aux integration variables
	VECTOR q_old;//variables in previous time-step
	VECTOR qvel_old;
	VECTOR qaccel_old;

	VECTOR _internalForces;
	VECTOR _externalForces;

	COO_MATRIX dampingMatrix;
	COO_MATRIX _stiffnessMatrix;
	COO_MATRIX systemMatrix;
	//these two store the damping parameters
	double dampingMassCoef;
	double dampingStiffnessCoef;

	double timestep;

	TetMesh *_tetMesh;
	Material *_material;
public:
	FullIntegrator(TetMesh *tetMesh)
	{ 
		_tetMesh = tetMesh; 
		_material = new MaterialLinearElasticity(1e7, 0.4);
		_internalForces = _material->getInternalForce();
		_stiffnessMatrix = _material->getStiffnessMatrix();
	}
	~FullIntegrator(){}

	virtual void timeStepSystem();
};
#endif