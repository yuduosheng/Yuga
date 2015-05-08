#include "StVKInternalForces.h"

StVKInternalForces::StVKInternalForces(TetMesh * volumetricMesh_, StVKTetABCD * precomputedABCDIntegrals_, bool addGravity_, double g_) : volumetricMesh(volumetricMesh_), precomputedIntegrals(precomputedABCDIntegrals_), gravityForce(NULL), addGravity(addGravity_), g(g_)
{
	numElements = volumetricMesh->total_tetrahedra;
	numVertices = volumetricMesh->total_points;
	numElementVertices = 4;

	lambdaLame = volumetricMesh->lambda;
	muLame = volumetricMesh->mu;


	buffer.resize(volumetricMesh->total_points);
	numElementVertices = volumetricMesh->total_points;
	InitGravity();
}

StVKInternalForces::~StVKInternalForces()
{
	gravityForce.clear();
	buffer.clear();
}

void StVKInternalForces::InitGravity()
{
	if (addGravity)
	{
		gravityForce.resize(volumetricMesh->total_points);
		computeGravity(&gravityForce, g);
	}
}
void StVKInternalForces::computeGravity(vector<dvec3>* gravityForce, double g, bool addForce=true) const
{
	if (!addForce)
		memset(gravityForce, 0, sizeof(dvec3) * numVertices);

	double invNumElementVertices = 1.0 / numElementVertices;

	for (int el = 0; el < numElements; el++)
	{
		double volume = volumetricMesh->tetrahedra[el].volume;
		double density = 1000;
		double mass = density * volume;
		for (int j = 0; j < numElementVertices; j++)
			gravityForce[volumetricMesh->tetrahedra[el].indices[j]][1] -= invNumElementVertices * mass * g; // gravity assumed to act in negative y-direction
	}
}
double StVKInternalForces::ComputeEnergy(vector<dvec3>* vertexDisplacements)
{
	return ComputeEnergyContribution(vertexDisplacements, 0, numElements);
}

double StVKInternalForces::ComputeEnergyContribution(vector<dvec3>* vertexDisplacements, int elementLow, int elementHigh)
{

	double energy = 0;

	// ---- linear
	memset(&buffer[0], 0, sizeof(dvec3) * numVertices);

	AddLinearTermsContribution(vertexDisplacements, &buffer, elementLow, elementHigh);
	for (int i = 0; i< 3 * numVertices; i++)
		energy += 0.5 * dot(buffer[i], (*vertexDisplacements)[i]);

	// ---- quadratic
	memset(&buffer[0], 0, sizeof(dvec3) * numVertices);
	AddQuadraticTermsContribution(vertexDisplacements, &buffer, elementLow, elementHigh);
	double oneThird = 1.0 / 3;
	for (int i = 0; i< 3 * numVertices; i++)
		energy += oneThird * dot(buffer[i], (*vertexDisplacements)[i]);

	// ---- cubic
	memset(&buffer[0], 0, sizeof(dvec3) * numVertices);
	AddCubicTermsContribution(vertexDisplacements, &buffer, elementLow, elementHigh);
	double oneQuarter = 1.0 / 4;
	for (int i = 0; i< 3 * numVertices; i++)
		energy += oneQuarter * dot(buffer[i], (*vertexDisplacements)[i]);

	return energy;
}

void StVKInternalForces::ComputeForces(vector<dvec3> * vertexDisplacements, vector<dvec3>* forces)
{
	//PerformanceCounter forceCounter;
	memset(forces, 0, sizeof(dvec3) * numVertices);

	AddLinearTermsContribution(vertexDisplacements, forces);
	AddQuadraticTermsContribution(vertexDisplacements, forces);
	AddCubicTermsContribution(vertexDisplacements, forces);

	if (addGravity)
	{
		int n = numVertices;
		for (int i = 0; i<3 * n; i++)
			(*forces)[i] -= gravityForce[i];
	}

	//forceCounter.StopCounter();
	//printf("Internal forces: %G\n", forceCounter.GetElapsedTime());
}

void StVKInternalForces::AddLinearTermsContribution(vector<dvec3>* vertexDisplacements, vector<dvec3>* forces, int elementLow, int elementHigh)
{
	if (elementLow < 0)
		elementLow = 0;
	if (elementHigh < 0)
		elementHigh = numElements;

	int * vertices = (int*)malloc(sizeof(int) * numElementVertices);

	void * elIter;
	precomputedIntegrals->AllocateElementIterator(&elIter);

	for (int el = elementLow; el < elementHigh; el++)
	{
		precomputedIntegrals->PrepareElement(el, elIter);
		for (int ver = 0; ver<numElementVertices; ver++)
			vertices[ver] = volumetricMesh->getVertexIndex(el, ver);

		double lambda = lambdaLame;
		double mu = muLame;

		for (int c = 0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
		{
			// linear terms
			for (int a = 0; a<numElementVertices; a++) // over all vertices
			{
				dvec3 qa(vertexDisplacements[3 * vertices[a] + 0],
					vertexDisplacements[3 * vertices[a] + 1],
					vertexDisplacements[3 * vertices[a] + 2]);

				dvec3 force = lambda * (precomputedIntegrals->A(elIter, c, a) * qa) +
					(mu * precomputedIntegrals->B(elIter, a, c)) * qa +
					mu * (precomputedIntegrals->A(elIter, a, c) * qa);

				(*forces)[vertices[c]] += force;

			}
		}
	}

	free(vertices);

	precomputedIntegrals->ReleaseElementIterator(elIter);
}

void StVKInternalForces::AddQuadraticTermsContribution(vector<dvec3>* vertexDisplacements, vector<dvec3>* forces, int elementLow, int elementHigh)
{
	if (elementLow < 0)
		elementLow = 0;
	if (elementHigh < 0)
		elementHigh = numElements;

	int * vertices = (int*)malloc(sizeof(int) * numElementVertices);

	void * elIter;
	precomputedIntegrals->AllocateElementIterator(&elIter);

	for (int el = elementLow; el < elementHigh; el++)
	{
		precomputedIntegrals->PrepareElement(el, elIter);
		for (int ver = 0; ver<numElementVertices; ver++)
			vertices[ver] = volumetricMesh->getVertexIndex(el, ver);

		double lambda = lambdaLame;
		double mu = muLame;

		for (int c = 0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
		{
			// quadratic terms
			for (int a = 0; a<numElementVertices; a++) // over all vertices
				for (int b = 0; b<numElementVertices; b++)
				{
					/*
					Vec3d force(0,0,0);
					Vec3d qa(vertexDisplacements[3*vertices[a]+0],
					vertexDisplacements[3*vertices[a]+1],
					vertexDisplacements[3*vertices[a]+2]);

					Vec3d qb(vertexDisplacements[3*vertices[b]+0],
					vertexDisplacements[3*vertices[b]+1],
					vertexDisplacements[3*vertices[b]+2]);

					double dotp = dot(qa,qb);

					force += 0.5 * lambda * dotp * precomputedIntegrals->C(el,c,a,b) +
					mu * dotp * precomputedIntegrals->C(el,a,b,c);

					Vec3d C = lambda * precomputedIntegrals->C(el,a,b,c) +
					mu * (precomputedIntegrals->C(el,c,a,b) + precomputedIntegrals->C(el,b,a,c));

					force += dot(C,qa) * qb;

					forces[3*vertices[c]+0] += force[0];
					forces[3*vertices[c]+1] += force[1];
					forces[3*vertices[c]+2] += force[2];
					*/

					dvec3 qa = (*vertexDisplacements)[vertices[a]];

					dvec3 qb = (*vertexDisplacements)[vertices[b]];

					double dotp = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2];

					dvec3 forceTerm1 = 0.5 * lambda * dotp * precomputedIntegrals->C(elIter, c, a, b) +
						mu * dotp * precomputedIntegrals->C(elIter, a, b, c);

					dvec3 C = lambda * precomputedIntegrals->C(elIter, a, b, c) +
						mu * (precomputedIntegrals->C(elIter, c, a, b) + precomputedIntegrals->C(elIter, b, a, c));

					double dotCqa = C[0] * qa[0] + C[1] * qa[1] + C[2] * qa[2];

					(*forces)[vertices[c]] +=  dvec3(forceTerm1[0] + dotCqa * qb[0],
												  forceTerm1[1] + dotCqa * qb[1],
												  forceTerm1[2] + dotCqa * qb[2]);
				}
		}
	}

	free(vertices);

	precomputedIntegrals->ReleaseElementIterator(elIter);
}

void StVKInternalForces::AddCubicTermsContribution(vector<dvec3>* vertexDisplacements, vector<dvec3>* forces, int elementLow, int elementHigh)
{
	if (elementLow < 0)
		elementLow = 0;
	if (elementHigh < 0)
		elementHigh = numElements;

	int * vertices = (int*)malloc(sizeof(int) * numElementVertices);

	void * elIter;
	precomputedIntegrals->AllocateElementIterator(&elIter);

	for (int el = elementLow; el < elementHigh; el++)
	{
		precomputedIntegrals->PrepareElement(el, elIter);

		for (int ver = 0; ver<numElementVertices; ver++)
			vertices[ver] = volumetricMesh->getVertexIndex(el, ver);

		double lambda = lambdaLame;
		double mu = muLame;

		for (int c = 0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
		{
			int vc = vertices[c];
			// cubic terms
			for (int a = 0; a<numElementVertices; a++) // over all vertices
			{
				int va = vertices[a];
				for (int b = 0; b<numElementVertices; b++)
				{
					int vb = vertices[b];
					for (int d = 0; d<numElementVertices; d++)
					{
						int vd = vertices[d];
						/*
						Vec3d qa(vertexDisplacements[3*va+0],
						vertexDisplacements[3*va+1],
						vertexDisplacements[3*va+2]);

						Vec3d qb(vertexDisplacements[3*vb+0],
						vertexDisplacements[3*vb+1],
						vertexDisplacements[3*vb+2]);

						Vec3d qd(vertexDisplacements[3*vd+0],
						vertexDisplacements[3*vd+1],
						vertexDisplacements[3*vd+2]);

						double dotp = dot(qa,qb);

						Vec3d force = 0.5 * lambda * dotp * precomputedIntegrals_->D(a,b,c,d) * qd +
						mu * dotp * precomputedIntegrals_->D(a,c,b,d) * qd;

						forces[3*vertices[c]+0] += force[0];
						forces[3*vertices[c]+1] += force[1];
						forces[3*vertices[c]+2] += force[2];
						*/
						dvec3 qa = (*vertexDisplacements)[va];
						dvec3 qb = (*vertexDisplacements)[vb];
						dvec3 qd = (*vertexDisplacements)[vd];

						double dotp = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2];
						double scalar = dotp * (0.5 * lambda * precomputedIntegrals->D(elIter, a, b, c, d) + mu * precomputedIntegrals->D(elIter, a, c, b, d));

						(*forces)[vc] += dvec3( scalar * qd[0],
												scalar * qd[1],
												scalar * qd[2]);
						
					}
				}
			}
		}
	}

	free(vertices);

	precomputedIntegrals->ReleaseElementIterator(elIter);
}
