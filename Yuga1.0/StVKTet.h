#ifndef _STVKTET_H_
#define _STVKTET_H_

#include "tetMesh.h"

class StvkTet
{
public:
	StvkTet();
	StvkTet(TetMesh * mesh);
	~StvkTet();

	void Precomputation(vector<glm::dmat3> *B);
	void ComputeElasticForces(vector<glm::dvec3> *f);
	void ComputeForceDifferentials(vector<glm::dvec3> detaX, vector<glm::dvec3>  *detaF);
private:
	float nu_ = 0.33f;			//Poisson ratio
	float E_ = 15000.0f;		//Young modulus
	double lambda;
	double mu;


	TetMesh * tetmesh;

	vector<glm::dvec3> deltaX;

	vector<glm::dmat3> B;
	vector<double> W;
	vector<glm::dvec3> f;
	vector<glm::dvec3> detaF;
};



#endif