#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
using namespace std;

typedef struct
{
	double volume;
	glm::dvec3 Phig[4]; // gradient of a basis function
} elementData;

void StVKSingleTetABCD(glm::dvec3 vtx[4], elementData * target)
{
	//double det = TetMesh::getTetDeterminant(&vtx[0], &vtx[1], &vtx[2], &vtx[3]);
	double det = 1.0f;
	target->volume = 0.33f;//fabs(det / 6);

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
					cout <<countJ<<" "<<countI<< " " <<ii << " " << jj << endl;
					countJ++;
				}
				countI++;
			}
			int sign = (((i + j) % 2) == 0) ? 1 : -1;
			target->Phig[i][j] = 1.0 * sign * glm::dot(glm::dvec3(1, 1, 1), glm::cross(columns[0], columns[1])) / det;
			
		}
}

int main()
{
	glm::dvec3 vtx[4];
	vtx[0] = glm::dvec3(1.0f, 0, 0);
	vtx[1] = glm::dvec3(0.0f, 1.0f, 0);
	vtx[2] = glm::dvec3(0.0f, 0, 1.0f);
	vtx[3] = glm::dvec3(0.0f, 0, 0);
	elementData el;
	StVKSingleTetABCD(vtx, &el);

	cout << el.Phig[0][0] << " " << el.Phig[0][1] << " " << el.Phig[0][1] << endl;
	cout << el.Phig[1][0] << " " << el.Phig[1][1] << " " << el.Phig[1][1] << endl;
	cout << el.Phig[2][0] << " " << el.Phig[2][1] << " " << el.Phig[2][1] << endl;
	cout << el.Phig[3][0] << " " << el.Phig[3][1] << " " << el.Phig[3][1] << endl;

	system("pause");
	return 0;
}