/*
This class is a container for a tetrahedral volumetric 3D mesh.
*/

#ifndef _TETMESH_H_
#define _TETMESH_H_

#include "SETTINGS.h"
#include "COO_MATRIX.h"

struct BoundaryTriangle
{
	UINT indices[3];
};

struct Tetrahedron {
	int indices[4];			//indices
	double volume;			//volume 
	//float plastic[6];		//plasticity values
	VEC3F e1, e2, e3;	//edges
	//glm::mat3 Re;			//Rotational warp of tetrahedron.
	//glm::mat3 Ke[4][4];		//Stiffness element matrix
	//glm::vec3 B[4];			//Jacobian of shapefunctions; B=SN =[d/dx  0     0 ][wn 0  0]
	//                                  [0    d/dy   0 ][0 wn  0]
	//									[0     0   d/dz][0  0 wn]
	//									[d/dy d/dx   0 ]
	//									[d/dz  0   d/dx]
	//									[0    d/dz d/dy]
};
class TetMesh
{
public:
	// loads a file of a "special" (not .veg) type
	// currently one such special format is supported:
	// specialFileType=0: 
	//   the ".ele" and ".node" format, used by TetGen, 
	//   "filename" is the basename, e.g., passing "mesh" will load the mesh from "mesh.ele" and "mesh.node" 
	// default material parameters will be used
	TetMesh(const char * filename);
	void ReadModelFromFile(const char * filename);
	// === misc queries ===
	void computeMassMatrix();

	VEC3F getElementCenter(int el) const;

	double GetTetraVolume(VEC3F e1, VEC3F e2, VEC3F e3) 
	{
		return e1.dot(e2.cross(e3)) / 6.0f;
	}

	static double getTetDeterminant(VEC3F * a, VEC3F * b, VEC3F * c, VEC3F * d);

	double getElementVolume(int el);

	void getElementInertiaTensor(int el, MATRIX3 & inertiaTensor);
	void computeElementMassMatrix(int element, double * massMatrix);
	//=====render=======
	void CalculateVN();
	void RenderModel();

	//====member queries====
	int GetVetexNumber(){ return total_points; }
	int GetTetNumber(){ return total_tetrahedra; }
	int getVertexIndex(int el, int ver);

protected:

	void AddBTriangle(int i0, int i1, int i2)
	{
		BoundaryTriangle t;
		t.indices[0] = i0;
		t.indices[1] = i1;
		t.indices[2] = i2;

		bTriangle.push_back(t);
	}
	void AddTetrahedron(int i0, int i1, int i2, int i3) {
		Tetrahedron t;

		t.indices[0] = i0;
		t.indices[1] = i1;
		t.indices[2] = i2;
		t.indices[3] = i3;

		t.e1 = oriCoordinates[i1] - oriCoordinates[i0];
		t.e2 = oriCoordinates[i2] - oriCoordinates[i0];
		t.e3 = oriCoordinates[i3] - oriCoordinates[i0];

		t.volume = GetTetraVolume(t.e1, t.e2, t.e3);

		tetrahedra.push_back(t);
	}

	//void computeElementMassMatrixHelper(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d, double * buffer);

protected:
	vector<VEC3F> oriCoordinates;
	vector<VEC3F> curCoordinates;

	vector<Tetrahedron> tetrahedra;
	vector<BoundaryTriangle> bTriangle;

	vector<VEC3F> normalOfTriangle;
	vector<VEC3F> normalOfVetex;
	vector<vector<int>> trianglesOfVertices;

	int total_points;
	int total_tetrahedra;
	int total_btriangle;

	GLuint                  meshVBuffer;
	GLuint                  meshNBuffer;
	GLuint                  indiceBuffer;
};

#endif

