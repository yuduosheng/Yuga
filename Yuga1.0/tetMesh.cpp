#include "tetMesh.h"

ofstream debug("debug.txt");

TetMesh::TetMesh(const char * filename)
{
	lambda = nu_ * E_ / ((1 + nu_)*(1 - 2 * nu_));
	mu = 0.5 * E_ / (1 + nu_);

	ReadModelFromFile(filename);
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

	oriCoordinates.resize(total_points);
	curCoordinates.resize(total_points);
	normalOfVetex.resize(total_points);
	trianglesOfVertices.resize(total_points);

	int noUse;
	nodeFile >> noUse >> noUse >> noUse;

	glm::vec3 nodePosition;

	for (int i = 0; i < total_points; ++i)
	{
		nodeFile >> noUse;
		nodeFile >> nodePosition.x >> nodePosition.y >> nodePosition.z;
		oriCoordinates[i] = nodePosition;
		curCoordinates[i] = oriCoordinates[i];
	}

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

void TetMesh::computeElementMassMatrix(int el, double * massMatrix) const
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
glm::dvec3 TetMesh::getElementCenter(int el) const
{
	glm::dvec3 pos(0, 0, 0);
	for (int i = 0; i < 4; ++i)
		pos += curCoordinates[tetrahedra[el].indices[i]];

	pos *= 0.25;

	return pos;
}

double TetMesh::getElementVolume(int el) const
{

}

void TetMesh::getElementInertiaTensor(int el, glm::dmat3 & inertiaTensor) const
{
	glm::dvec3 a = curCoordinates[tetrahedra[el].indices[0]];
	glm::dvec3 b = curCoordinates[tetrahedra[el].indices[1]];
	glm::dvec3 c = curCoordinates[tetrahedra[el].indices[2]];
	glm::dvec3 d = curCoordinates[tetrahedra[el].indices[3]];

	glm::dvec3 center = getElementCenter(el);
	a -= center;
	b -= center;
	c -= center;
	d -= center;

	double absdetJ = 6.0 * getElementVolume(el);

	double x1 = a[0], x2 = b[0], x3 = c[0], x4 = d[0];
	double y1 = a[1], y2 = b[1], y3 = c[1], y4 = d[1];
	double z1 = a[2], z2 = b[2], z3 = c[2], z4 = d[2];

	double A = absdetJ * (y1*y1 + y1*y2 + y2*y2 + y1*y3 + y2*y3 + y3*y3 + y1*y4 + y2*y4 + y3*y4 + y4*y4 + z1*z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

	double B = absdetJ * (x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3 + x1*x4 + x2*x4 + x3*x4 + x4*x4 + z1*z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

	double C = absdetJ * (x1*x1 + x1*x2 + x2*x2 + x1*x3 + x2*x3 + x3*x3 + x1*x4 + x2*x4 + x3*x4 + x4*x4 + y1*y1 + y1 * y2 + y2 * y2 + y1 * y3 + y2 * y3 + y3 * y3 + y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4) / 60.0;

	double Ap = absdetJ * (2 * y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2 + 2 * y2*z2 + y3*z2 + y4*z2 + y1*z3 + y2*z3 + 2 * y3*z3 + y4*z3 + y1*z4 + y2*z4 + y3*z4 + 2 * y4*z4) / 120.0;

	double Bp = absdetJ * (2 * x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2 + 2 * x2*z2 + x3*z2 + x4*z2 + x1*z3 + x2*z3 + 2 * x3*z3 + x4*z3 + x1*z4 + x2*z4 + x3*z4 + 2 * x4*z4) / 120.0;

	double Cp = absdetJ * (2 * x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2 + 2 * x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + 2 * x3*y3 + x4*y3 + x1*y4 + x2*y4 + x3*y4 + 2 * x4*y4) / 120.0;

	inertiaTensor = glm::dmat3(A, -Bp, -Cp, -Bp, B, -Ap, -Cp, -Ap, C);
}

bool TetMesh::containsVertex(int el, glm::dvec3 pos) const // true if given element contain given position, false otherwise
{
	double weights[4];
	computeBarycentricWeights(el, pos, weights);

	// all weights must be non-negative
	return ((weights[0] >= 0) && (weights[1] >= 0) && (weights[2] >= 0) && (weights[3] >= 0));
}

void TetMesh::computeBarycentricWeights(int el, glm::dvec3 pos, double * weights) const
{
	/*
	|x1 y1 z1 1|
	D0 = |x2 y2 z2 1|
	|x3 y3 z3 1|
	|x4 y4 z4 1|

	|x  y  z  1|
	D1 = |x2 y2 z2 1|
	|x3 y3 z3 1|
	|x4 y4 z4 1|

	|x1 y1 z1 1|
	D2 = |x  y  z  1|
	|x3 y3 z3 1|
	|x4 y4 z4 1|

	|x1 y1 z1 1|
	D3 = |x2 y2 z2 1|
	|x  y  z  1|
	|x4 y4 z4 1|

	|x1 y1 z1 1|
	D4 = |x2 y2 z2 1|
	|x3 y3 z3 1|
	|x  y  z  1|

	wi = Di / D0
	*/

	glm::dvec3 vtx[4];
	for (int i = 0; i<4; i++)
		vtx[i] = curCoordinates[tetrahedra[el].indices[i]];

	double D[5];
	D[0] = getTetDeterminant(&vtx[0], &vtx[1], &vtx[2], &vtx[3]);

	for (int i = 1; i <= 4; i++)
	{
		glm::dvec3 buf[4];
		for (int j = 0; j<4; j++)
			buf[j] = vtx[j];
		buf[i - 1] = pos;
		D[i] = getTetDeterminant(&buf[0], &buf[1], &buf[2], &buf[3]);
		weights[i - 1] = D[i] / D[0];
	}
}

double TetMesh::getTetDeterminant(glm::dvec3 * a, glm::dvec3 * b, glm::dvec3 * c, glm::dvec3 * d)
{
	// computes the determinant of the 4x4 matrix
	// [ a 1 ]
	// [ b 1 ]
	// [ c 1 ]
	// [ d 1 ]
	glm::dvec4 aa = glm::dvec4(a->x, b->x, c->x, d->x);
	glm::dvec4 bb = glm::dvec4(a->y, b->y, c->z, d->y);
	glm::dvec4 cc = glm::dvec4(a->z, b->y, c->z, d->z);
	glm::dvec4 dd = glm::dvec4(1, 1, 1, 1);

	glm::mat4 mat = glm::mat4(aa, bb, cc, dd);

	return glm::determinant(mat);
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
		glm::vec3 e01, e02;
		e01 = curCoordinates[n1] - curCoordinates[n0];
		e02 = curCoordinates[n2] - curCoordinates[n0];
		normalOfTriangle[i] = glm::normalize(glm::cross(e02, e01));
	}

	for (int i = 0; i < total_points; ++i)
	{
		glm::vec3 normal(0);
		for (int j = 0; j < trianglesOfVertices[i].size(); ++j)
		{
			normal += normalOfTriangle[trianglesOfVertices[i][j]];
		}
		normalOfVetex[i] = glm::normalize(normal);
	}
}

void TetMesh::RenderModel()
{
	CalculateVN();
	
	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		curCoordinates.size() * sizeof(glm::vec3),
		&(curCoordinates[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	glGenBuffers(1, &meshNBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		normalOfVetex.size() * sizeof(glm::vec3),
		&(normalOfVetex[0]),
		GL_STATIC_DRAW);

	// 2rst attribute buffer : normal
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);


	glGenBuffers(1, &indiceBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bTriangle.size() * sizeof(BoundaryTriangle), &(bTriangle[0]), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);

	// Draw the triangles !
	glDrawElements(
		GL_TRIANGLES,      // mode
		bTriangle.size() * 3,    // count
		GL_UNSIGNED_SHORT,   // type
		(void*)0           // element array buffer offset
		);
	//draw vertices normal
	/*
	glBegin(GL_LINES);
	for (int i = 0; i < total_points; ++i)
	{
		glVertex3f(curCoordinates[i].x, curCoordinates[i].y, curCoordinates[i].z);
		glVertex3f(normalOfVetex[i].x + curCoordinates[i].x, normalOfVetex[i].y + curCoordinates[i].y, normalOfVetex[i].z + curCoordinates[i].z);
	}
	glEnd();*/
	
}

/*
void TetMesh::interpolateGradient(int element, const double * U, int numFields, glm::vec3 pos, double * grad) const
{
	computeGradient(element, U, numFields, grad);
}

void TetMesh::computeGradient(int el, const double * U, int numFields, double * grad) const
{
	// grad is 9 x numFields
	// grad is constant inside a tet
	glm::vec3 vtx[4];
	for (int i = 0; i<4; i++)
		vtx[i] = curCoordinates[tetrahedra[el].indices[i]];

	// form M =
	// [b - a]
	// [c - a]
	// [d - a]

	Mat3d M(vtx[1] - vtx[0], vtx[2] - vtx[0], vtx[3] - vtx[0]);
	Mat3d MInv = inv(M);
	//printf("M=\n");
	//M.print();

	for (int field = 0; field<numFields; field++)
	{
		// form rhs =
		// [U1 - U0]
		// [U2 - U0]
		// [U3 - U0]
		const double * u[4];
		for (int i = 0; i<4; i++)
			u[i] = &U[3 * numVertices * field + 3 * getVertexIndex(element, i)];

		glm::vec3 rows[3];
		for (int i = 0; i<3; i++)
			rows[i] = glm::vec3(u[i + 1]) - glm::vec3(u[0]);

		Mat3d rhs(rows);
		//printf("rhs=\n");
		//rhs.print();

		Mat3d gradM = trans(MInv * rhs);
		gradM.convertToArray(&grad[9 * field]);
	}
}

void TetMesh::orient()
{
	for (int el = 0; el<total_tetrahedra; el++)
	{
		double det = dot(*(getVertex(el, 1)) - *(getVertex(el, 0)), cross(*(getVertex(el, 2)) - *(getVertex(el, 0)), *(getVertex(el, 3)) - *(getVertex(el, 0))));

		if (det < 0)
		{
			// reverse tet
			int * elementVertices = elements[el];
			// swap 2 and 3
			int buff = elementVertices[2];
			elementVertices[2] = elementVertices[3];
			elementVertices[3] = buff;
		}
	}
}
*/

