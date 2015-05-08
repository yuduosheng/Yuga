#include "corotationalLinearFEM.h"

CorotationalLinearFEM::CorotationalLinearFEM(TetMesh * tetMesh_)
{
	tetMesh = tetMesh_;
	numVertices = tetMesh->total_points;


	F0.resize(numVertices);
	mass.resize(numVertices);
	A_row.resize(total_points);
	K_row.resize(total_points);


	// store the undeformed positions
	int numElements = tetMesh->total_tetrahedra;

	// set Lame's coefficients
	float d15 = E_ / (1.0f + nu_) / (1.0f - 2 * nu_);
	float d16 = (1.0f - nu_) * d15;
	float d17 = nu_ * d15;
	float d18 = E_ / 2 / (1.0f + nu_);

	for (int i = 0; i < tetMesh->total_tetrahedra; ++i)
	{
		glm::vec3 e10 = tetMesh->tetrahedra[i].e1;
		glm::vec3 e20 = tetMesh->tetrahedra[i].e2;
		glm::vec3 e30 = tetMesh->tetrahedra[i].e3;
		//Eq. 10.32
		glm::mat3 E = glm::mat3(e10.x, e10.y, e10.z,
			e20.x, e20.y, e20.z,
			e30.x, e30.y, e30.z);
		float detE = glm::determinant(E);
		float invDetE = 1.0f / detE;

		//Eq. 10.40 (a) & Eq. 10.42 (a)
		//Shape function derivatives wrt x,y,z
		// d/dx N0
		float invE10 = (e20.z*e30.y - e20.y*e30.z)*invDetE;
		float invE20 = (e10.y*e30.z - e10.z*e30.y)*invDetE;
		float invE30 = (e10.z*e20.y - e10.y*e20.z)*invDetE;
		float invE00 = -invE10 - invE20 - invE30;

		//Eq. 10.40 (b) & Eq. 10.42 (b)
		// d/dy N0
		float invE11 = (e20.x*e30.z - e20.z*e30.x)*invDetE;
		float invE21 = (e10.z*e30.x - e10.x*e30.z)*invDetE;
		float invE31 = (e10.x*e20.z - e10.z*e20.x)*invDetE;
		float invE01 = -invE11 - invE21 - invE31;

		//Eq. 10.40 (c) & Eq. 10.42 (c)
		// d/dz N0
		float invE12 = (e20.y*e30.x - e20.x*e30.y)*invDetE;
		float invE22 = (e10.x*e30.y - e10.y*e30.x)*invDetE;
		float invE32 = (e10.y*e20.x - e10.x*e20.y)*invDetE;
		float invE02 = -invE12 - invE22 - invE32;
/*
		tetMesh->tetrahedra[i].B[0] = glm::vec3(invE00, invE01, invE02);
		tetMesh->tetrahedra[i].B[1] = glm::vec3(invE10, invE11, invE12);
		tetMesh->tetrahedra[i].B[2] = glm::vec3(invE20, invE21, invE22);
		tetMesh->tetrahedra[i].B[3] = glm::vec3(invE30, invE31, invE32);

		for (int i = 0; i<4; i++) {
			for (int j = 0; j<4; j++) {
				glm::mat3& Ke = tetMesh->tetrahedra[i].Ke[i][j];
				float d19 = tetMesh->tetrahedra[i].B[i].x;
				float d20 = tetMesh->tetrahedra[i].B[i].y;
				float d21 = tetMesh->tetrahedra[i].B[i].z;
				float d22 = tetMesh->tetrahedra[i].B[j].x;
				float d23 = tetMesh->tetrahedra[i].B[j].y;
				float d24 = tetMesh->tetrahedra[i].B[j].z;
				Ke[0][0] = d16 * d19 * d22 + d18 * (d20 * d23 + d21 * d24);
				Ke[0][1] = d17 * d19 * d23 + d18 * (d20 * d22);
				Ke[0][2] = d17 * d19 * d24 + d18 * (d21 * d22);

				Ke[1][0] = d17 * d20 * d22 + d18 * (d19 * d23);
				Ke[1][1] = d16 * d20 * d23 + d18 * (d19 * d22 + d21 * d24);
				Ke[1][2] = d17 * d20 * d24 + d18 * (d21 * d23);

				Ke[2][0] = d17 * d21 * d22 + d18 * (d19 * d24);
				Ke[2][1] = d17 * d21 * d23 + d18 * (d20 * d24);
				Ke[2][2] = d16 * d21 * d24 + d18 * (d20 * d23 + d19 * d22);

				Ke *= tetMesh->tetrahedra[i].volume;
			}
		}*/
	}
}

CorotationalLinearFEM::~CorotationalLinearFEM()
{

}

void CorotationalLinearFEM::RecalcMassMatrix() {
	//This is a lumped mass matrix
	//Based on Eq. 10.106 and pseudocode in Fig. 10.9 on page 358
	for (size_t i = 0; i<tetMesh->total_points; i++) {
		mass[i] = 1.0f / tetMesh->total_points;
	}

	for (int i = 0; i<tetMesh->total_tetrahedra; i++) {
		float m = (density*tetMesh->tetrahedra[i].volume)* 0.25f;
		mass[tetMesh->tetrahedra[i].indices[0]] += m;
		mass[tetMesh->tetrahedra[i].indices[1]] += m;
		mass[tetMesh->tetrahedra[i].indices[2]] += m;
		mass[tetMesh->tetrahedra[i].indices[3]] += m;
	}
}
