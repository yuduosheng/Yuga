#ifndef MATERIAL_H
#define MATERIAL_H

#include "SETTINGS.h"

class Material
{
protected:
	Real _Young;
	Real _Poisson;
	Real _mu;
	Real _lambda;

	VECTOR _solutionX;
	VECTOR _InternalForce;


	Material(Real Young, Real Poisson)
	{
		_Young = Young;
		_Poisson = Poisson;
		_mu = 0.5 * _Young / (1 + _Poisson);
		_lambda = _Young * _Poisson / (1 + _Poisson) * (1 - _Poisson - _Poisson);
	}
	~Material() {}

	virtual MATRIX3 firstPiolaKirchhoff(const MATRIX3 &F) = 0;
	virtual VECTOR& computeInternalForce() = 0;

};

#endif