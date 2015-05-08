#include "StVKTet.h"

StvkTet::StvkTet(TetMesh * mesh)
{
	tetmesh = mesh;
	lambda = nu_ * E_ / ((1 + nu_)*(1-2*nu_));
	mu = 0.5 * E_ / (1 + nu_);
	B.resize(tetmesh->total_tetrahedra);
	f.resize(tetmesh->total_points);
	deltaX.resize(tetmesh->total_points);

}
StvkTet::~StvkTet()
{
	B.clear();
	f.clear();
	deltaX.clear();
}

void StvkTet::Precomputation(vector<glm::dmat3> *B)
{
	glm::dmat3 Dm;
	for (int tn = 0; tn < tetmesh->total_tetrahedra; ++tn)
	{
		glm::dvec3 c0 = glm::dvec3(tetmesh->tetrahedra[tn].e1.x, tetmesh->tetrahedra[tn].e2.x, tetmesh->tetrahedra[tn].e3.x);
		glm::dvec3 c1 = glm::dvec3(tetmesh->tetrahedra[tn].e1.y, tetmesh->tetrahedra[tn].e2.y, tetmesh->tetrahedra[tn].e3.y);
		glm::dvec3 c2 = glm::dvec3(tetmesh->tetrahedra[tn].e1.z, tetmesh->tetrahedra[tn].e2.z, tetmesh->tetrahedra[tn].e3.z);
		Dm[0] = c0;
		Dm[1] = c1;
		Dm[2] = c2;

		(*B)[tn] = glm::inverse(Dm);
	}
}
void StvkTet::ComputeElasticForces(vector<glm::dvec3> *f)
{
	memset(f, 0, sizeof(*f));
	for (int tn = 0; tn < tetmesh->total_tetrahedra; ++tn)
	{
		glm::dmat3 Ds;
		int i, j, k, l;
		i = tetmesh->tetrahedra[tn].indices[1];
		j = tetmesh->tetrahedra[tn].indices[2];
		k = tetmesh->tetrahedra[tn].indices[3];
		l = tetmesh->tetrahedra[tn].indices[0];
		Ds[0] = glm::dvec3(tetmesh->curCoordinates[i].x - tetmesh->curCoordinates[l].x,
			tetmesh->curCoordinates[j].x - tetmesh->curCoordinates[l].x, 
			tetmesh->curCoordinates[k].x - tetmesh->curCoordinates[l].x);

		Ds[1] = glm::dvec3(tetmesh->curCoordinates[i].y - tetmesh->curCoordinates[l].y,
			tetmesh->curCoordinates[j].y - tetmesh->curCoordinates[l].y,
			tetmesh->curCoordinates[k].y - tetmesh->curCoordinates[l].y);

		Ds[2] = glm::dvec3(tetmesh->curCoordinates[i].z - tetmesh->curCoordinates[l].z,
			tetmesh->curCoordinates[j].z - tetmesh->curCoordinates[l].z,
			tetmesh->curCoordinates[k].z - tetmesh->curCoordinates[l].z);

		glm::dmat3 F = Ds * B[tn];

		glm::dmat3 E = 0.5 * (glm::transpose(F)*F) - 0.5 * glm::dmat3(1.0);
		glm::dmat3 P = 2 * mu * F * E;
		double trE = E[0][0] + E[1][1] + E[2][2];

		P += lambda * trE * F *glm::dmat3(1);


		glm::dmat3 H = -1.0 * tetmesh->tetrahedra[tn].volume * P * (glm::transpose(B[tn]));

		(*f)[i] += glm::dvec3(H[0].x, H[1].x, H[2].x);
		(*f)[j] += glm::dvec3(H[0].y, H[1].y, H[2].y);
		(*f)[k] += glm::dvec3(H[0].z, H[1].z, H[2].z);
		(*f)[l] -= glm::dvec3(H[0].x, H[1].x, H[2].x) + glm::dvec3(H[0].y, H[1].y, H[2].y) + glm::dvec3(H[0].z, H[1].z, H[2].z);

	}
}

void StvkTet::ComputeForceDifferentials(vector<glm::dvec3> detaX, vector<glm::dvec3>  *detaf)
{
	memset(detaf, 0, sizeof(*detaf));

	for (int tn = 0; tn < tetmesh->total_tetrahedra; ++tn)
	{
		glm::dmat3 Ds;
		glm::dmat3 detaDs;
		int i, j, k, l;
		i = tetmesh->tetrahedra[tn].indices[1];
		j = tetmesh->tetrahedra[tn].indices[2];
		k = tetmesh->tetrahedra[tn].indices[3];
		l = tetmesh->tetrahedra[tn].indices[0];
		Ds[0] = glm::dvec3(tetmesh->curCoordinates[i].x - tetmesh->curCoordinates[l].x,
			tetmesh->curCoordinates[j].x - tetmesh->curCoordinates[l].x,
			tetmesh->curCoordinates[k].x - tetmesh->curCoordinates[l].x);

		Ds[1] = glm::dvec3(tetmesh->curCoordinates[i].y - tetmesh->curCoordinates[l].y,
			tetmesh->curCoordinates[j].y - tetmesh->curCoordinates[l].y,
			tetmesh->curCoordinates[k].y - tetmesh->curCoordinates[l].y);

		Ds[2] = glm::dvec3(tetmesh->curCoordinates[i].z - tetmesh->curCoordinates[l].z,
			tetmesh->curCoordinates[j].z - tetmesh->curCoordinates[l].z,
			tetmesh->curCoordinates[k].z - tetmesh->curCoordinates[l].z);

		detaDs[0] = glm::dvec3(detaX[i].x - detaX[l].x, detaX[j].x - detaX[l].x, detaX[k].x - detaX[l].x);
		detaDs[1] = glm::dvec3(detaX[i].y - detaX[l].y, detaX[j].y - detaX[l].y, detaX[k].y - detaX[l].y);
		detaDs[2] = glm::dvec3(detaX[i].z - detaX[l].z, detaX[j].z - detaX[l].z, detaX[k].z - detaX[l].z);

		glm::dmat3 F = Ds * B[tn];
		glm::dmat3 detaF = detaDs * B[tn];
		glm::dmat3 E = 0.5 * (glm::transpose(F)*F) - 0.5 * glm::dmat3(1.0);
		glm::dmat3 detaE = 0.5 * glm::transpose(detaF) * F + 0.5 * glm::transpose(F) * detaF;
		double trE = E[0][0] + E[1][1] + E[2][2];
		double trdetaE = detaE[0][0] + detaE[1][1] + detaE[2][2];
		glm::dmat3 detaP = 2 * mu * (detaF * E + F * detaE);
		detaP += lambda * trE * detaF * glm::dmat3(1);
		detaP += lambda * trdetaE * F * glm::dmat3(1);

		glm::dmat3 detaH = -1.0 * tetmesh->tetrahedra[tn].volume * detaP * (glm::transpose(B[tn]));

		(*detaf)[i] += glm::dvec3(detaH[0].x, detaH[1].x, detaH[2].x);
		(*detaf)[j] += glm::dvec3(detaH[0].y, detaH[1].y, detaH[2].y);
		(*detaf)[k] += glm::dvec3(detaH[0].z, detaH[1].z, detaH[2].z);
		(*detaf)[l] -= glm::dvec3(detaH[0].x, detaH[1].x, detaH[2].x) + glm::dvec3(detaH[0].y, detaH[1].y, detaH[2].y) + glm::dvec3(detaH[0].z, detaH[1].z, detaH[2].z);
	}
}