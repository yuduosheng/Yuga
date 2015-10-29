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
#include "CONJUGATE_GRADIENT.h"
extern ofstream debug;
CONJUGATE_GRADIENT::CONJUGATE_GRADIENT(bool pre)
{
	_maxIteration = 1000;//SIMPLE_PARSER::getInt("cg max iterations", 10);
	_eps = 0.01;// SIMPLE_PARSER::getFloat("cg eps", 0.0001);

	if (pre)
		_preconditioner = new JACOBI_PRECONDITIONER;

}

bool CONJUGATE_GRADIENT::solve(COO_MATRIX &A, const VECTOR& b, VECTOR& x)
{
	_residual.resize(A.rows());
	_direction.resize(A.rows());
	_q.resize(A.rows());
	_s.resize(A.rows());

	if (b.size() == 0)
		return true;

	if (_preconditioner == NULL)
		return solveCG(A, b, x);
	else
	{
		_preconditioner->setDiag(&A.getDiag());
		return solvePCG(A, b, x);
	}
		
}

bool CONJUGATE_GRADIENT::solveCG(COO_MATRIX &A, const VECTOR& b, VECTOR& x)
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
		//cout << "Conjugate Gradient solver converged in " << i << " iterations" << endl;
	}
	else
		cout << "Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;
	return converged;
}


bool CONJUGATE_GRADIENT::solvePCG(COO_MATRIX &A, const VECTOR& b, VECTOR& x)
{
	//VECTOR diag = _preconditioner->getDiag();
	//debug << diag << endl;
	//debug << "=================-------------------" << endl;
	//vector<TRIPLET>& matrix = A.matrix();
	//cout << matrix.size() << " " << A.rows() * A.rows() << endl;
	//for (int i = 0; i < matrix.size(); i++)
	//{
	//	debug << matrix[i].row() << "  " << matrix[i].col() << "  " << matrix[i].value() << endl;
	//}
	//cout << endl;

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
		//cout << "Conjugate Gradient solver converged in " << i << " iterations" << endl;
	}
	else
		cout << "Preconditioned Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;

	return converged;
}