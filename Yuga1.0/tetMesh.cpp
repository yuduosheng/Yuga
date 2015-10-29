#include "tetMesh.h"

extern ofstream debug;

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
	//debug << bTriangle.size() << endl;
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
		Tetrahedron& t = getTet(i);
		Dm[0] = oriCoordinates[t.indices[1]] - oriCoordinates[t.indices[0]];
		Dm[1] = oriCoordinates[t.indices[2]] - oriCoordinates[t.indices[0]];
		Dm[2] = oriCoordinates[t.indices[3]] - oriCoordinates[t.indices[0]];

		_InverDm[i].setZero();
		_InverDm[i] << Dm[0][0], Dm[1][0], Dm[2][0],
					   Dm[0][1], Dm[1][1], Dm[2][1],
					   Dm[0][2], Dm[1][2], Dm[2][2];

		_InverDm[i] = _InverDm[i].inverse().eval();
		
	}

	//compute PFPu
	for (int i = 0; i < total_tetrahedra; ++i)
	{
		computePFPu(tetrahedra[i]._PFPu, _InverDm[i]);
		//debug << tetrahedra[i]._PFPu << endl;
		//debug << "==========================================" << endl;
	}

	computeMassMatrix();
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

void TetMesh::setCurPosition(VECTOR & q)
{
	for (int i = 0; i < total_points; i++)
    {
    	curCoordinates[i] += q.segment<3>(3 * i);
    }
}
void TetMesh::getDisplacement(VECTOR &q)
{
	for (int i = 0; i < total_points; i++)
    {
    	q.segment<3>(3 * i) = curCoordinates[i] - oriCoordinates[i];
    }
}
int TetMesh::getClosestVertex(VEC3F pos)
{
	// linear scan
	double closestDist = DBL_MAX;
	int closestVertex = -1;

	for (int i = 0; i < total_points; i++)
	{
		VEC3F & vertexPosition = curCoordinates[i];
		double dist = (pos - vertexPosition).norm();
		if (dist < closestDist)
		{
			closestDist = dist;
			closestVertex = i;
		}
	}

	return closestVertex;
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

void TetMesh::computePFPu(MATRIX& _PFPu, MATRIX3& matInv)//p(F)/p(x)
{
	const double m = matInv(0, 0);
	const double n = matInv(0, 1);
	const double o = matInv(0, 2);
	const double p = matInv(1, 0);
	const double q = matInv(1, 1);
	const double r = matInv(1, 2);
	const double s = matInv(2, 0);
	const double t = matInv(2, 1);
	const double u = matInv(2, 2);

	const double t1 = -m - p - s;
	const double t2 = -n - q - t;
	const double t3 = -o - r - u;

	_PFPu.resize(9, 12);
	_PFPu.setZero();
	_PFPu(0, 0) = t1;                  //t1  0   0   m    0    0   p   0   0   s   0   0
	_PFPu(0, 1) = 0.0;                 //0   t1  0   0    m    0   0   p   0   0   s   0
	_PFPu(0, 2) = 0.0;                 //0   0   t1  0    0    m   0   0   p   0   0   s
	_PFPu(0, 3) = m;                   //t2  0   0   n    0    0   q   0   0   t   0   0   
	_PFPu(0, 4) = 0.0;                 //0   t2  0   0    n    0   0   q   0   0   t   0
	_PFPu(0, 5) = 0.0;                 //0   0   t2  0    0    n   0   0   q   0   0   t
	_PFPu(0, 6) = p;                   //t3  0   0   o    0    0   r   0   0   u   0   0
	_PFPu(0, 7) = 0.0;                 //0   t3  0   0    o    0   0   r   0   0   u   0
	_PFPu(0, 8) = 0.0;                 //0   0   t3  0    0    o   0   0   r   0   0   u
	_PFPu(0, 9) = s;
	_PFPu(0, 10) = 0.0;
	_PFPu(0, 11) = 0.0;
	_PFPu(1, 0) = 0.0;
	_PFPu(1, 1) = t1;
	_PFPu(1, 2) = 0.0;
	_PFPu(1, 3) = 0.0;
	_PFPu(1, 4) = m;
	_PFPu(1, 5) = 0.0;
	_PFPu(1, 6) = 0.0;
	_PFPu(1, 7) = p;
	_PFPu(1, 8) = 0.0;
	_PFPu(1, 9) = 0.0;
	_PFPu(1, 10) = s;
	_PFPu(1, 11) = 0.0;
	_PFPu(2, 0) = 0.0;
	_PFPu(2, 1) = 0.0;
	_PFPu(2, 2) = t1;
	_PFPu(2, 3) = 0.0;
	_PFPu(2, 4) = 0.0;
	_PFPu(2, 5) = m;
	_PFPu(2, 6) = 0.0;
	_PFPu(2, 7) = 0.0;
	_PFPu(2, 8) = p;
	_PFPu(2, 9) = 0.0;
	_PFPu(2, 10) = 0.0;
	_PFPu(2, 11) = s;
	_PFPu(3, 0) = t2;
	_PFPu(3, 1) = 0.0;
	_PFPu(3, 2) = 0.0;
	_PFPu(3, 3) = n;
	_PFPu(3, 4) = 0.0;
	_PFPu(3, 5) = 0.0;
	_PFPu(3, 6) = q;
	_PFPu(3, 7) = 0.0;
	_PFPu(3, 8) = 0.0;
	_PFPu(3, 9) = t;
	_PFPu(3, 10) = 0.0;
	_PFPu(3, 11) = 0.0;
	_PFPu(4, 0) = 0.0;
	_PFPu(4, 1) = t2;
	_PFPu(4, 2) = 0.0;
	_PFPu(4, 3) = 0.0;
	_PFPu(4, 4) = n;
	_PFPu(4, 5) = 0.0;
	_PFPu(4, 6) = 0.0;
	_PFPu(4, 7) = q;
	_PFPu(4, 8) = 0.0;
	_PFPu(4, 9) = 0.0;
	_PFPu(4, 10) = t;
	_PFPu(4, 11) = 0.0;
	_PFPu(5, 0) = 0.0;
	_PFPu(5, 1) = 0.0;
	_PFPu(5, 2) = t2;
	_PFPu(5, 3) = 0.0;
	_PFPu(5, 4) = 0.0;
	_PFPu(5, 5) = n;
	_PFPu(5, 6) = 0.0;
	_PFPu(5, 7) = 0.0;
	_PFPu(5, 8) = q;
	_PFPu(5, 9) = 0.0;
	_PFPu(5, 10) = 0.0;
	_PFPu(5, 11) = t;
	_PFPu(6, 0) = t3;
	_PFPu(6, 1) = 0.0;
	_PFPu(6, 2) = 0.0;
	_PFPu(6, 3) = o;
	_PFPu(6, 4) = 0.0;
	_PFPu(6, 5) = 0.0;
	_PFPu(6, 6) = r;
	_PFPu(6, 7) = 0.0;
	_PFPu(6, 8) = 0.0;
	_PFPu(6, 9) = u;
	_PFPu(6, 10) = 0.0;
	_PFPu(6, 11) = 0.0;
	_PFPu(7, 0) = 0.0;
	_PFPu(7, 1) = t3;
	_PFPu(7, 2) = 0.0;
	_PFPu(7, 3) = 0.0;
	_PFPu(7, 4) = o;
	_PFPu(7, 5) = 0.0;
	_PFPu(7, 6) = 0.0;
	_PFPu(7, 7) = r;
	_PFPu(7, 8) = 0.0;
	_PFPu(7, 9) = 0.0;
	_PFPu(7, 10) = u;
	_PFPu(7, 11) = 0.0;
	_PFPu(8, 0) = 0.0;
	_PFPu(8, 1) = 0.0;
	_PFPu(8, 2) = t3;
	_PFPu(8, 3) = 0.0;
	_PFPu(8, 4) = 0.0;
	_PFPu(8, 5) = o;
	_PFPu(8, 6) = 0.0;
	_PFPu(8, 7) = 0.0;
	_PFPu(8, 8) = r;
	_PFPu(8, 9) = 0.0;
	_PFPu(8, 10) = 0.0;
	_PFPu(8, 11) = u;
}