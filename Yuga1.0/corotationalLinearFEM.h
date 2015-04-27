#ifndef _COROTATIONALLINEARFEM_H_
#define _COROTATIONALLINEARFEM_H_

/*
Corotational linear FEM deformable model.

This class implements the deformable model described in the following paper:

M. Mueller, M. Gross: Interactive Virtual Materials.
In Proc. of Graphics Interface 2004 (2004), pp. 239–246.

In [Mueller 2004], the tangent stiffness matrix is approximate (warp=1).
This class can also compute the exact tangent stiffness matrix (warp=2).
The implementation is described in:
J. Barbic: Exact Corotational Linear FEM Stiffness Matrix, Technical Report, USC, 2012

It is also possible to turn warping off (warp=0). This gives fast linear FEM dynamics,
but large deformations are not well-represented.
*/

#include "tetMesh.h"
#include <map>

class CorotationalLinearFEM
{
public:

	// initializes corotational linear FEM
	// input: tetMesh
	CorotationalLinearFEM(TetMesh * tetMesh);
	virtual ~CorotationalLinearFEM();

	// computes the internal forces and (warped) stiffness matrix for the entire mesh
	// vertex displacements (input) and internal forces (output) must be (pre-allocated) vectors of length 3 * numVertices
	// the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
	// i.e., the computed internal forces are *negatives* of the actual physical internal forces acting on the material
	// warp:
	//   0: no warping (linear FEM)
	//   1: stiffness warping (corotational linear FEM with approximate stiffness matrix) [Mueller 2004]
	void RecalcMassMatrix();
	void InitializePlastic();
	void ClearStiffnessAssembly();
	void InitializePlastic();



	inline TetMesh * GetTetMesh() { return tetMesh; }

protected:
	float nu_ = 0.33f;			//Poisson ratio
	float E_ = 15000.0f;		//Young modulus
	float density = 1000.0f;
	int numVertices;
	TetMesh * tetMesh;

	vector<float> mass;			//Mass matrix
	vector<glm::vec3> F0;

	typedef map<int, glm::mat3> matrix_map;
	typedef matrix_map::iterator matrix_iterator;
	vector<matrix_map> K_row;
	vector<matrix_map> A_row;
};

#endif

