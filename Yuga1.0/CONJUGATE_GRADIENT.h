/******************************************************************************
 * This file is part of Cubica++.
 
 * Cubica++ is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * Cubica++ is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details.

 * You should have received a copy of the GNU General Public License along
 * with Cubica++. If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "SETTINGS.h"
#include <iostream>

using namespace::std;

class JACOBI_PRECONDITIONER
{
public:
	JACOBI_PRECONDITIONER(VECTOR& diagonal) :
		_diag(diagonal)
	{

	}
	~JACOBI_PRECONDITIONER() {};

	void solve(const VECTOR& b, VECTOR& x)
	{
		x = (b.array() / _diag.array()).matrix();
	}

private:
	const VECTOR& _diag;
};


class CONJUGATE_GRADIENT
{
public:
	CONJUGATE_GRADIENT(COO_MATRIX& A_);
  ~CONJUGATE_GRADIENT();

  int& maxIteration() { return _maxIteration; };
  Real& eps()   { return _eps; };

  bool solve(const VECTOR& b, VECTOR& x);
  bool solveCG(const VECTOR& b, VECTOR& x);
  bool solvePCG(const VECTOR& b, VECTOR& x);

  vector<int>& numberOfIterations() { return _numberOfIterations; };

protected:
  JACOBI_PRECONDITIONER * _preconditioner;

  int _maxIteration;
  Real _eps;

  int _totalFrames;
  int _totalIterations;
  vector<int> _numberOfIterations;

  COO_MATRIX A;

  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  VECTOR _s;
};
/******************************************************************************
* This file is part of Cubica++.

* Cubica++ is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free
* Software Foundation, either version 3 of the License, or (at your option)
* any later version.

* Cubica++ is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.

* You should have received a copy of the GNU General Public License along
* with Cubica++. If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

CONJUGATE_GRADIENT::CONJUGATE_GRADIENT(COO_MATRIX& A_) :
_totalFrames(0),
_totalIterations(0)
{
	_maxIteration = 1000;//SIMPLE_PARSER::getInt("cg max iterations", 10);
	_eps = 0.01;// SIMPLE_PARSER::getFloat("cg eps", 0.0001);
	A = A_;
	_preconditioner = new JACOBI_PRECONDITIONER(A.getDiag());

}

CONJUGATE_GRADIENT::~CONJUGATE_GRADIENT()
{

}

bool CONJUGATE_GRADIENT::solve(const VECTOR& b, VECTOR& x)
{
	if (b.size() == 0)
		return true;

	if (_preconditioner == NULL)
		return solveCG(b, x);
	else
		return solvePCG(b, x);
}

bool CONJUGATE_GRADIENT::solveCG(const VECTOR& b, VECTOR& x)
{
	x.conservativeResize(b.size());
	x.setZero();

	//_integrator->getMatVecMult(x, _q);
	A.multiplyVector(x, _q);

	_direction = _residual = b - _q;

	Real deltaNew = _residual.squaredNorm();

	Real delta0 = deltaNew;

	bool converged = false;
	int i = 0;
	for (; i < _maxIteration && !converged; i++) {

		//_integrator->getMatVecMult(_direction, _q);
		A.multiplyVector(_direction, _q);

		Real alpha = _direction.dot(_q);
		if (fabs(alpha) > 0.0)
			alpha = deltaNew / alpha;

		x += alpha * _direction;

		if (i > 0 && i % 25 == 0) {
			//_integrator->getMatVecMult(x, _residual);
			A.multiplyVector(x, _residual);
			_residual *= -1;
			_residual += b;
		}
		else
			_residual -= alpha * _q;

		Real deltaOld = deltaNew;

		deltaNew = _residual.squaredNorm();

		Real beta = deltaNew / deltaOld;

		_direction *= beta;
		_direction += _residual;

		if (deltaNew <= _eps * delta0)
			converged = true;
	}

	if (converged) {
		// cout << "Conjugate Gradient solver converged in " << i << " iterations" << endl;
		_totalIterations++;
		_totalFrames += i;
		_numberOfIterations.push_back(_totalFrames);
	}
	else
		cout << "Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;
	return converged;
}


bool CONJUGATE_GRADIENT::solvePCG(const VECTOR& b, VECTOR& x)
{
	x.conservativeResize(b.size());
	x.setZero();

	//_integrator->getMatVecMult(x, _q);
	A.multiplyVector(x, _q);
	_residual = b - _q;

	_preconditioner->solve(_residual, _direction);

	Real deltaNew = _residual.dot(_direction);

	Real delta0 = deltaNew;

	bool converged = false;
	int i = 0;
	for (; i < _maxIteration && !converged; i++) {
		//_integrator->getMatVecMult(_direction, _q);
		A.multiplyVector(_direction, _q);

		Real alpha = _direction.dot(_q);
		if (fabs(alpha) > 0.0)
			alpha = deltaNew / alpha;

		x += alpha * _direction;

		if (i > 0 && i % 25 == 0) {
			//_integrator->getMatVecMult(x, _residual);
			A.multiplyVector(x, _residual);
			_residual *= -1;
			_residual += b;
		}
		else
			_residual -= alpha * _q;

		_preconditioner->solve(_residual, _s);

		Real deltaOld = deltaNew;

		deltaNew = _residual.dot(_s);


		Real beta = deltaNew / deltaOld;

		_direction *= beta;
		_direction += _s;

		if (deltaNew <= _eps * delta0)
			converged = true;
	}

	if (converged) {
		_totalIterations++;
		_totalFrames += i;
		_numberOfIterations.push_back(i);
	}
	else
		cout << "Preconditioned Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;

	return converged;
}

#endif
