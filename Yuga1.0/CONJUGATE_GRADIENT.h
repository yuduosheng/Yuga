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
#include "COO_MATRIX.h"
#include <iostream>

using namespace::std;

class JACOBI_PRECONDITIONER
{
public:
	JACOBI_PRECONDITIONER(){ _diag = NULL; }
	JACOBI_PRECONDITIONER(VECTOR& diagonal) :
		_diag(&diagonal)
	{

	}
	~JACOBI_PRECONDITIONER() {};
	void setDiag(VECTOR* diag){ _diag = diag; }
	VECTOR& getDiag(){ return *_diag; };
	void solve(const VECTOR& b, VECTOR& x)
	{
		x = (b.array() / (*_diag).array()).matrix();
	}

private:
	VECTOR* _diag;
};


class CONJUGATE_GRADIENT
{
public:
	CONJUGATE_GRADIENT(bool pre);
	~CONJUGATE_GRADIENT()
	{
		if (_preconditioner)
			delete _preconditioner;
	};

  int& maxIteration() { return _maxIteration; };
  Real& eps()   { return _eps; };

  bool solve(COO_MATRIX &A, const VECTOR& b, VECTOR& x);
  bool solveCG(COO_MATRIX &A, const VECTOR& b, VECTOR& x);
  bool solvePCG(COO_MATRIX &A, const VECTOR& b, VECTOR& x);

protected:
  JACOBI_PRECONDITIONER * _preconditioner;

  int _maxIteration;
  Real _eps;

  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  VECTOR _s;
};

#endif
