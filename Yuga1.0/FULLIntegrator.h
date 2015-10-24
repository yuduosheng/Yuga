#ifndef FULLINTEGRATOR_H
#define FULLINTEGRATOR_H
#include "TetMesh.h"
#include "Material.h"
#include "COO_MATRIX.h"
class FullIntegrator
{
private:
	VECTOR displacement;
	TetMesh &_tetMesh;
public:
	FullIntegrator(TetMesh &tetMesh)
	{ 
		_tetMesh = tetMesh; 
	}
	~FullIntegrator(){}
};
#endif