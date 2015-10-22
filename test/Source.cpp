#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>


#include <iostream>

using namespace std;

typedef Eigen::Vector3d VEC3F;
typedef Eigen::Vector2d VEC2F;
typedef Eigen::VectorXd VECTOR;
typedef Eigen::ArrayXd ARRAY;
typedef Eigen::MatrixXd MATRIX;
typedef double Real;
typedef Eigen::Matrix3d MATRIX3;
typedef Eigen::Matrix4d MATRIX4;
typedef Eigen::Matrix<double, 3, 4, Eigen::ColMajor> MATRIX3x4;
typedef Eigen::Matrix<double, 9, 9, Eigen::RowMajor> MATRIX9;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Quaternion<double> QUATERNION;
typedef Eigen::Triplet<double> TRIPLET;
typedef Eigen::Vector3i VEC3I;
typedef unsigned int UINT;

#include "COO_MATRIX.h"

int main()
{
	//MATRIX4 mat;
	//mat <<  2, 1, 1, 1,
	//		1, 2, 1, 1,
	//		1, 1, 2, 1,
	//		1, 1, 1, 2;
	MATRIX3 mat3;
	//cout << mat3 << endl;
	//cout << mat << endl;
	//COO_MATRIX mat;
	//mat.add(1.1, 0, 0);
	//mat.add(1.2, 0, 1);
	//mat.add(1.3, 1, 0);
	//mat.add(1.3, 1, 1);
	//
	//mat.add(1.1, 1, 0);
	//mat.add(1.1, 1, 1);
	//
	//mat.add(2, 2, 2);
	//mat.PrintEntry();

	mat3 << 3, 0, 0,
		0, 4, 0,
		0, 0, 3;
	cout << mat3.inverse().eval() << endl;
	system("pause");
	return 0;
}