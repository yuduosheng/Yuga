#ifndef _STVKINTERNALFORCES_H_
#define _STVKINTERNALFORCES_H_

#include "tetMesh.h"
#include "StVKTetABCD.h"

class StVKInternalForces
{
public:

	// before creating this class, you must first create the volumetric mesh, and the precomputed integrals
	// you can use the StVKElementABCDLoader.h header file to create the "precomputedABCDIntegrals"
	StVKInternalForces(TetMesh * tetMesh, StVKTetABCD * precomputedABCDIntegrals, bool addGravity = false, double g = 9.81);
	virtual ~StVKInternalForces();

	// both vertex displacements and internal forces refer to the vertices of the simulation mesh
	// they must be (pre-allocated) vectors of length 3 * numVertices
	// the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
	// i.e., the computed internal forces are negatives of the actual physical internal forces acting on the material
	virtual void ComputeForces(vector<dvec3>* vertexDisplacements, vector<dvec3>* internalForces);

	// enables or disables the gravity (note: you can also set this in the constructor; use this routine to turn the gravity on/off during the simulation)
	void SetGravity(bool addGravity) { this->addGravity = addGravity; InitGravity(); } // if addGravity is enabled, ComputeForces will subtract the gravity force from the internal forces (note: subtraction, not addition, is used because the internal forces are returned with the sign as described in the f_int(x) comment above)

	void computeGravity(vector<dvec3> * gravityForce, double g, bool addForce) const;

	virtual double ComputeEnergy(vector<dvec3>* vertexDisplacements); // get the nonlinear elastic strain energy

	inline TetMesh * GetVolumetricMesh() { return volumetricMesh; }
	inline StVKTetABCD * GetPrecomputedIntegrals() { return precomputedIntegrals; }



	// === advanced routines below === 
	double ComputeEnergyContribution(vector<dvec3> * vertexDisplacements, int elementLow, int elementHigh); // compute the contribution to strain energy due to the specified elements; needs a buffer for internal calculations; you can pass NULL (and then an internal buffer will be used), or pass your own buffer (useful with multi-threading)
	void AddLinearTermsContribution(vector<dvec3> * vertexDisplacements, vector<dvec3>* forces, int elementLow = -1, int elementHigh = -1);
	void AddQuadraticTermsContribution(vector<dvec3>* vertexDisplacements, vector<dvec3>* forces, int elementLow = -1, int elementHigh = -1);
	void AddCubicTermsContribution(vector<dvec3>* vertexDisplacements, vector<dvec3>* forces, int elementLow = -1, int elementHigh = -1);

protected:
	TetMesh * volumetricMesh;
	StVKTetABCD * precomputedIntegrals;

	vector<dvec3> gravityForce;
	bool addGravity;
	double g;
	void InitGravity(); // aux function

	vector<dvec3> buffer;
	int numElementVertices;
	int numVertices;
	int numElements;

	double lambdaLame;
	double muLame;
};

#endif