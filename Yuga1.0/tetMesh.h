/*
This class is a container for a tetrahedral volumetric 3D mesh.
*/

#ifndef _TETMESH_H_
#define _TETMESH_H_

#include "App.h"


struct BoundaryTriangle
{
	GLushort indices[3];
};

struct Tetrahedron {
	int indices[4];			//indices
	float volume;			//volume 
	float plastic[6];		//plasticity values
	glm::vec3 e1, e2, e3;	//edges
	glm::mat3 Re;			//Rotational warp of tetrahedron.
	glm::mat3 Ke[4][4];		//Stiffness element matrix
	glm::vec3 B[4];			//Jacobian of shapefunctions; B=SN =[d/dx  0     0 ][wn 0  0]
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
	glm::dvec3 getElementCenter(int el) const;
	static double getTetVolume(glm::vec3 * a, glm::vec3 * b, glm::vec3 * c, glm::vec3 * d);
	static double getTetDeterminant(glm::vec3 * a, glm::vec3 * b, glm::vec3 * c, glm::vec3 * d);
	double getElementVolume(int el) const;
	void getElementInertiaTensor(int el, glm::dmat3 & inertiaTensor) const;
	void computeElementMassMatrix(int element, double * massMatrix) const;

	bool containsVertex(int element, glm::vec3 pos) const; // true if given element contain given position, false otherwise

	// === interpolation ===

	void computeBarycentricWeights(int el, glm::vec3 pos, double * weights) const;
	//void computeGradient(int element, const double * U, int numFields, double * grad) const; // for tet meshes, gradient is constant inside each tet, hence no need to specify position
	//void interpolateGradient(int element, const double * U, int numFields, glm::vec3 pos, double * grad) const; // conforms to the virtual function in the base class, "pos" does not affect the computation

	// === misc ===

	//void orient(); // orients the tets (re-orders vertices within each tet), so that each tet has positive orientation: ((v1 - v0) x (v2 - v0)) dot (v3 - v0) >= 0
	

	//=====render=======
	void CalculateVN();
	void RenderModel();
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

		tetrahedra.push_back(t);
	}

	//void computeElementMassMatrixHelper(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d, double * buffer);

private:
	vector<glm::vec3> oriCoordinates;
	vector<glm::vec3> curCoordinates;

	vector<Tetrahedron> tetrahedra;
	vector<BoundaryTriangle> bTriangle;

	vector<glm::vec3> normalOfTriangle;
	vector<glm::vec3> normalOfVetex;
	vector<vector<int>> trianglesOfVertices;

	int total_points = 0;
	int total_tetrahedra = 0;
	int total_btriangle = 0;

	GLuint                  meshVBuffer;
	GLuint                  meshNBuffer;
	GLuint                  indiceBuffer;
};

#endif
