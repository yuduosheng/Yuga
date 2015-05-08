#ifndef _STVKTETABCD_H_
#define _STVKTETABCD_H_

#include "tetMesh.h"

/*
This class stores the St.Venant-Kirchhoff A,B,C,D coefficients for a tetrahedral element.
This is the low-memory version (the version that we use most often).
See also StVKInternalForces.h .
*/

class StVKTetABCD
{
public:

	// computes the ABCD coefficients 
	StVKTetABCD(TetMesh * tetMesh);

	glm::dmat3 A(void * elementIterator, int i, int j);
	double B(void * elementIterator, int i, int j);
	glm::dvec3 C(void * elementIterator, int i, int j, int k);
	double D(void * elementIterator, int i, int j, int k, int l);

	typedef struct
	{
		double volume;
		glm::dvec3 Phig[4]; // gradient of a basis function
	} elementData;

	typedef struct
	{
		elementData * elementPointer;
		double dots[4][4];
	} ElementCache;

	void AllocateElementIterator(void ** elementIterator);
	void ReleaseElementIterator(void * elementIterator);
	void PrepareElement(int el, void * elementIterator); // must call each time before accessing an element

	virtual ~StVKTetABCD();


	inline dmat3 tensorProduct(dvec3 & vecA, dvec3 & vecB)
	{
		dmat3 result(vecA[0] * vecB[0], vecA[0] * vecB[1], vecA[0] * vecB[2],
			vecA[1] * vecB[0], vecA[1] * vecB[1], vecA[1] * vecB[2],
			vecA[2] * vecB[0], vecA[2] * vecB[1], vecA[2] * vecB[2]);

		return result;
	}

protected:

	elementData * elementsData;

	// creates the elementData structure for a tet
	void StVKSingleTetABCD(glm::dvec3 vertices[4], elementData * target);
};

#endif
