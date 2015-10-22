#ifndef MATERIALLINEARELASTICITY_H
#define MATERIALLINEARELASTICITY_H
#include "Material.h"
#include "TetMesh.h"

class MaterialLinearElasticity : public Material
{
private:
	TetMesh *_tetMesh;

public:
	MaterialLinearElasticity(Real Young, Real Poisson) : Material(Young, Poisson)
	{
		_InternalForce.resize(_tetMesh->GetVetexNumber() * 3);
	};
	~MaterialLinearElasticity(){};

	MATRIX3 firstPiolaKirchhoff(const MATRIX3 &F)
	{
		return _mu * (F + F.transpose() - MATRIX3::Identity() - MATRIX3::Identity()) + _lambda * (F - MATRIX3::Identity()).trace() * MATRIX3::Identity();
	}

	VECTOR& computeInternalForce()
	{
		_InternalForce.setZero();

		vector<VEC3F>& curVer = _tetMesh->getCurVertexBuffer();
		for (int i = 0; i < _tetMesh->GetTetNumber; ++i)
		{
			Tetrahedron& t = _tetMesh->getTet(i);
			MATRIX3 F;
			VEC3F rowDs[3];
			rowDs[0] = curVer[t.indices[1]] - curVer[t.indices[0]];
			rowDs[1] = curVer[t.indices[2]] - curVer[t.indices[0]];
			rowDs[2] = curVer[t.indices[3]] - curVer[t.indices[0]];

			F.setZero();
			F << rowDs[0][0], rowDs[1][0], rowDs[2][0],
				rowDs[0][1], rowDs[1][1], rowDs[2][1],
				rowDs[0][2], rowDs[1][2], rowDs[2][2];

			MATRIX3& Bm = _tetMesh->getInverseDm(i);

			F = F * Bm;
			MATRIX3 P;
			P = firstPiolaKirchhoff(F);
			P = -t.volume * P * Bm.transpose();

			VEC3F forces[4];
			forces[1][0] = P(0, 0);    forces[2][0] = P(0, 1);	   forces[3][0] = P(0, 2);
			forces[1][1] = P(1, 0);    forces[2][1] = P(1, 1);	   forces[3][1] = P(1, 2);
			forces[1][2] = P(2, 0);	   forces[2][2] = P(2, 1);	   forces[3][2] = P(2, 2);

			forces[0] = -forces[1] - forces[2] - forces[3];

			for (int j = 0; j < 4; ++j)
			{
				_InternalForce[t.indices[i] * 3] += forces[i][0];
				_InternalForce[t.indices[i] * 3 + 1] += forces[i][1];
				_InternalForce[t.indices[i] * 3 + 2] += forces[i][2];
			}
		}
	}
};
#endif