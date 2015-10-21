#ifndef MATERIALLINEARELASTICITY_H
#define MATERIALLINEARELASTICITY_H
#include "Material.h"

class MaterialLinearElasticity : public Material
{
public:
	MaterialLinearElasticity(Real Young, Real Poisson) : Material(Young, Poisson)
	{};
	~MaterialLinearElasticity(){};

	MATRIX3 firstPiolaKirchhoff(const MATRIX3 &F)
	{
		return _mu * (F + F.transpose() - MATRIX3::Identity() - MATRIX3::Identity()) + _lambda * (F - MATRIX3::Identity()).trace() * MATRIX3::Identity();
	}
};
#endif