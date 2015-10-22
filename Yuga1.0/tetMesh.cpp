#include "tetMesh.h"

ofstream debug("debug.txt");

TetMesh::TetMesh(const char * filename)
{
	ReadModelFromFile(filename);
	init();
}

void TetMesh::ReadModelFromFile(const char *filename)
{
	string fileName = filename;
	string suffix = ".node";
	string suffix2 = ".ele";
	string suffix3 = ".face";

	fileName += suffix;
	ifstream nodeFile(fileName.c_str(), ios::in);
	if (!nodeFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	nodeFile >> total_points;

	normalOfVetex.resize(total_points);
	trianglesOfVertices.resize(total_points);

	int noUse;
	nodeFile >> noUse >> noUse >> noUse;

	VEC3F nodePosition;

	for (int i = 0; i < total_points; ++i)
	{
		nodeFile >> noUse;
		nodeFile >> nodePosition.x() >> nodePosition.y() >> nodePosition.z();
		oriCoordinates.push_back(nodePosition);
		curCoordinates.push_back(nodePosition);
		//debug << curCoordinates[i].x() <<" "<< curCoordinates[i].y() <<" "<< curCoordinates[i].z() << endl;
	}
	//debug << curCoordinates.size() << endl;
	nodeFile.close();

	fileName = filename;
	fileName += suffix2;

	ifstream eleFile(fileName.c_str(), ios::in);
	if (!eleFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	eleFile >> total_tetrahedra;
	eleFile >> noUse >> noUse;

	for (int i = 0; i < total_tetrahedra; ++i)
	{
		int p0, p1, p2, p3;
		eleFile >> noUse;
		eleFile >> p0 >> p1 >> p2 >> p3;
		AddTetrahedron(p0, p1, p2, p3);
	}
	eleFile.close();


	fileName = filename;
	fileName += suffix3;


	ifstream faceFile(fileName.c_str(), ios::in);
	if (!faceFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}
	faceFile >> total_btriangle;
	faceFile >> noUse;

	normalOfTriangle.resize(total_btriangle);

	for (int i = 0; i < total_btriangle; ++i)
	{
		int p0, p1, p2;
		faceFile >> noUse;
		faceFile >> p0 >> p1 >> p2;

		trianglesOfVertices[p0].push_back(i);
		trianglesOfVertices[p1].push_back(i);
		trianglesOfVertices[p2].push_back(i);

		AddBTriangle(p0, p1, p2);
	}
	debug << bTriangle.size() << endl;
	//for (int i = 0; i < bTriangle.size(); i++)
	//{
	//	debug << bTriangle[i].indices[0] <<"  "<< bTriangle[i].indices[1] << "  " <<bTriangle[i].indices[2] << endl;
	//}
	/*
	for (int i = 0; i < total_points; ++i)
	{
	for (int j = 0; j < trianglesOfVertices[i].size(); ++j)
	{
	debug << trianglesOfVertices[i][j]<<" ";
	}
	debug << endl;
	}
	*/

	faceFile.close();
}

void TetMesh::init()
{
	//compute tet InverDm
	VEC3F Dm[3];
	_InverDm.resize(total_tetrahedra);
	for (int i = 0; i < total_tetrahedra; ++i)
	{
		Dm[0] = oriCoordinates[getVertexIndex(i, 1)] - oriCoordinates[getVertexIndex(i, 0)];
		Dm[1] = oriCoordinates[getVertexIndex(i, 2)] - oriCoordinates[getVertexIndex(i, 0)];
		Dm[2] = oriCoordinates[getVertexIndex(i, 3)] - oriCoordinates[getVertexIndex(i, 0)];

		_InverDm[i].setZero();
		_InverDm[i] << Dm[0][0], Dm[1][0], Dm[2][0],
					   Dm[0][1], Dm[1][1], Dm[2][1],
					   Dm[0][2], Dm[1][2], Dm[2][2];

		_InverDm[i] = _InverDm[i].inverse();
	}
	
}

void TetMesh::computeElementMassMatrix(int el, double * massMatrix)
{
	/*
	Consistent mass matrix of a tetrahedron =

	[ 2  1  1  1  ]
	[ 1  2  1  1  ]
	mass / 20 * [ 1  1  2  1  ]
	[ 1  1  1  2  ]

	Note: mass = density * volume. Other than via the mass, the
	consistent mass matrix does not depend on the shape of the tetrahedron.
	(This can be seen after a long algebraic derivation; see:
	Singiresu S. Rao: The finite element method in engineering, 2004)
	*/

	const double mtx[16] = { 2, 1, 1, 1,
		1, 2, 1, 1,
		1, 1, 2, 1,
		1, 1, 1, 2 };

	double density = 1000;//getElementDensity(el);
	double factor = density * getElementVolume(el) / 20;

	for (int i = 0; i<16; i++)
		massMatrix[i] = factor * mtx[i];

	// lumped mass
	/*
	double mass = element(el)->density() * getElementVolume(el);
	massMatrix[0] = massMatrix[5] = massMatrix[10] = massMatrix[15] = mass / 4.0;
	*/
}
VEC3F TetMesh::getElementCenter(int el) const
{
	VEC3F pos(0, 0, 0);
	for (int i = 0; i < 4; ++i)
		pos += curCoordinates[tetrahedra[el].indices[i]];

	pos *= 0.25;

	return pos;
}

double TetMesh::getElementVolume(int el)
{
	VEC3F e1, e2, e3;	//edges
	int index[4];
	index[0] = tetrahedra[el].indices[0];
	index[1] = tetrahedra[el].indices[1];
	index[2] = tetrahedra[el].indices[2];
	index[3] = tetrahedra[el].indices[3];

	e1 = curCoordinates[index[1]] - curCoordinates[index[0]];
	e2 = curCoordinates[index[2]] - curCoordinates[index[0]];
	e3 = curCoordinates[index[3]] - curCoordinates[index[0]];
	
	return GetTetraVolume(e1, e2, e3);
}

double TetMesh::getTetDeterminant(VEC3F * a, VEC3F * b, VEC3F * c, VEC3F * d)
{
	// computes the determinant of the 4x4 matrix
	// [ a 1 ]
	// [ b 1 ]
	// [ c 1 ]
	// [ d 1 ]
	MATRIX4 mat;
	mat(0, 0) = a->x(); mat(0, 1) = a->y(); mat(0, 2) = a->z(); mat(0, 3) = 1.0f;
	mat(1, 0) = b->x(); mat(1, 1) = b->y(); mat(1, 2) = b->z(); mat(1, 3) = 1.0f;
	mat(2, 0) = c->x(); mat(2, 1) = c->y(); mat(2, 2) = c->z(); mat(2, 3) = 1.0f;
	mat(3, 0) = d->x(); mat(3, 1) = d->y(); mat(3, 2) = d->z(); mat(3, 3) = 1.0f;

	return mat.determinant();
}
int TetMesh::getVertexIndex(int el, int ver)
{
	return tetrahedra[el].indices[ver];
}

void TetMesh::CalculateVN()
{
	for (int i = 0; i < bTriangle.size(); ++i)
	{
		int n0 = bTriangle[i].indices[0];
		int n1 = bTriangle[i].indices[1];
		int n2 = bTriangle[i].indices[2];
		VEC3F e01, e02;
		e01 = curCoordinates[n1] - curCoordinates[n0];
		e02 = curCoordinates[n2] - curCoordinates[n0];
		//normalOfTriangle[i] = glm::normalize(glm::cross(e02, e01));
		normalOfTriangle[i] = e01.cross(e02).normalized();
	}

	for (int i = 0; i < total_points; ++i)
	{
		VEC3F normal(0, 0, 0);
		for (int j = 0; j < trianglesOfVertices[i].size(); ++j)
		{
			normal += normalOfTriangle[trianglesOfVertices[i][j]];
		}
		//normalOfVetex[i] = glm::normalize(normal);
		normalOfVetex[i] = normal.normalized();
	}
}

void TetMesh::RenderModel()
{
	CalculateVN();

	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		curCoordinates.size() * sizeof(VEC3F),
		&(curCoordinates[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

	glGenBuffers(1, &meshNBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		normalOfVetex.size() * sizeof(VEC3F),
		&(normalOfVetex[0]),
		GL_STATIC_DRAW);

	// 2rst attribute buffer : normal
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, NULL);


	glGenBuffers(1, &indiceBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bTriangle.size() * sizeof(BoundaryTriangle), &(bTriangle[0]), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);

	// Draw the triangles !
	glDrawElements(
		GL_TRIANGLES,      // mode
		bTriangle.size() * 3,    // count
		GL_UNSIGNED_INT,   // type
		(void*)0           // element array buffer offset
		);
	//draw vertices normal
	/*
	glBegin(GL_LINES);
	for (int i = 0; i < total_points; ++i)
	{
	glVertex3f(curCoordinates[i].x(), curCoordinates[i].y(), curCoordinates[i].z());
	glVertex3f(normalOfVetex[i].x() + curCoordinates[i].x(), normalOfVetex[i].y() + curCoordinates[i].y(), normalOfVetex[i].z() + curCoordinates[i].z());
	//glVertex3f(curCoordinates[i].x() + normalOfVetex[i].x(), curCoordinates[i].y() + normalOfVetex[i].y(), curCoordinates[i].z() + normalOfVetex[i].z());
	}
	glEnd();
	*/
}

void TetMesh::computeMassMatrix()
{
	//lumped mass
	double density = 1000.0;

	_massMatrix.resize(3 * total_points, 3 * total_points);

	for (int i = 0; i < total_tetrahedra; ++i)
	{
		double mass = density * getElementVolume(i) / 4.0f;
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 3; k++)
			{
				_massMatrix.add(mass, getVertexIndex(i, j) * 3 + k, getVertexIndex(i, j) * 3 + k);
			}	
	}
}