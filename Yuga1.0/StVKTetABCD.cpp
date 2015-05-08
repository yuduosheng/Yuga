#include "StVKTetABCD.h"

StVKTetABCD::StVKTetABCD(TetMesh * tetMesh)
{
	int numElements = tetMesh->total_tetrahedra;
	int totalCoefficientSize = (sizeof(elementData) * numElements);

	elementsData = (elementData*)malloc(sizeof(elementData) * numElements);

	for (int el = 0; el < numElements; el++)
	{
		glm::dvec3 vertices[4];
		for (int i = 0; i < 4; i++)
			vertices[i] = tetMesh->oriCoordinates[tetMesh->tetrahedra[el].indices[i]];
		StVKSingleTetABCD(vertices, &elementsData[el]);
	}

	printf("Total tet ABCD coefficient size: %G Mb.\n",
		1.0 * totalCoefficientSize / 1024 / 1024);
}

StVKTetABCD::~StVKTetABCD()
{
	free(elementsData);
}

void StVKTetABCD::StVKSingleTetABCD(glm::dvec3 vtx[4], elementData * target)
{
	double det = TetMesh::getTetDeterminant(&vtx[0], &vtx[1], &vtx[2], &vtx[3]);
	target->volume = fabs(det / 6);

	for (int i = 0; i<4; i++)
		for (int j = 0; j<3; j++)
		{
			glm::dvec3 columns[2];
			int countI = 0;
			for (int ii = 0; ii<4; ii++)
			{
				if (ii == i)
					continue;
				int countJ = 0;
				for (int jj = 0; jj<3; jj++)
				{
					if (jj == j)
						continue;

					columns[countJ][countI] = vtx[ii][jj];
					countJ++;
				}
				int sign = (((i + j) % 2) == 0) ? 1 : -1;
				target->Phig[i][j] = 1.0 * sign * glm::dot(glm::dvec3(1, 1, 1), glm::cross(columns[0], columns[1])) / det;
				countI++;
			}
		}
}

void StVKTetABCD::AllocateElementIterator(void ** elementIterator)
{
	ElementCache * cache = new ElementCache();
	*elementIterator = cache;
}

void StVKTetABCD::ReleaseElementIterator(void * elementIterator)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	delete(cache);
}

void StVKTetABCD::PrepareElement(int el, void * elementIterator)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	cache->elementPointer = &elementsData[el];
	for (int i = 0; i<4; i++)
		for (int j = 0; j<4; j++)
			(cache->dots)[i][j] = dot(cache->elementPointer->Phig[i], cache->elementPointer->Phig[j]);
}

dmat3 StVKTetABCD::A(void * elementIterator, int i, int j)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	return cache->elementPointer->volume * tensorProduct(cache->elementPointer->Phig[i], cache->elementPointer->Phig[j]);
}

double StVKTetABCD::B(void * elementIterator, int i, int j)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	return cache->elementPointer->volume * cache->dots[i][j];
}

glm::dvec3 StVKTetABCD::C(void * elementIterator, int i, int j, int k)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	return cache->elementPointer->volume * cache->dots[j][k] * cache->elementPointer->Phig[i];
}

double StVKTetABCD::D(void * elementIterator, int i, int j, int k, int l)
{
	ElementCache * cache = (ElementCache *)elementIterator;
	return cache->elementPointer->volume * cache->dots[i][j] * cache->dots[k][l];
}

