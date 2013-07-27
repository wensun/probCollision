#ifndef __MATRIX_H__
#define __MATRIX_H__

//#define _USE_MATH_DEFINES

#include <cmath>
#include <float.h>
#include <iostream>
#include <assert.h>
#include <limits>

template <size_t _numRows, size_t _numColumns = 1> class Matrix {

private:
	double _elems[_numRows * _numColumns];
	//double *_elems;

public:
	/*Matrix() { _elems = new double[_numRows * _numColumns]; }
	Matrix( const Matrix<_numRows, _numColumns>& q ) {
	_elems = new double[_numRows * _numColumns];
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
	_elems[i] = q._elems[i]; 
	}
	}

	~Matrix() { delete[] _elems; }

	Matrix<_numRows, _numColumns>& operator = (const Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
	_elems[i] = q._elems[i]; 
	}
	return (*this);
	}*/

	// Retrieval
	inline size_t numRows() const { 
		return _numRows; 
	}
	inline size_t numColumns() const { 
		return _numColumns; 
	}

	// Subscript operator
	inline double& operator () (size_t row, size_t column) {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}
	inline double  operator () (size_t row, size_t column) const {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}

	inline double& operator [] (size_t elt) {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}
	inline double  operator [] (size_t elt) const {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}

	// Reset to zeros
	inline void reset() { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] = double(0);
		}
	}

	// Submatrix
	template <size_t numRows, size_t numColumns>
	inline Matrix<numRows, numColumns> subMatrix(size_t row, size_t column) const {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		Matrix<numRows, numColumns> m;
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				m(i, j) = (*this)(row + i, column + j);
			}
		}
		return m;
	}

	// Insert
	template <size_t numRows, size_t numColumns>
	inline void insert(size_t row, size_t column, const Matrix<numRows, numColumns>& q) {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				(*this)(row + i, column + j) = q(i, j);
			}
		}
	}

	// Unary minus
	inline Matrix<_numRows, _numColumns> operator-() const {
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = -_elems[i];
		}
		return m;
	}

	// Unary plus
	inline const Matrix<_numRows, _numColumns>& operator+() const { 
		return *this; 
	}

	// Equality
	inline bool operator==(const Matrix<_numRows, _numColumns>& q) const { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
				return false;
			}
		}
		return true;
	}

	// Inequality
	inline bool operator!=(const Matrix<_numRows, _numColumns>& q) const { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
				return true;
			}
		}
		return false;
	}

	// Matrix addition
	inline Matrix<_numRows, _numColumns> operator+(const Matrix<_numRows, _numColumns>& q) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] + q._elems[i];
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator+=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] += q._elems[i];
		}
		return *this; 
	}

	// Matrix subtraction
	inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& q) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] - q._elems[i];
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator-=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] -= q._elems[i];
		}
		return *this; 
	}

	// Scalar multiplication
	inline Matrix<_numRows, _numColumns> operator*(double a) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] * a;
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator*=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline Matrix<_numRows, _numColumns> operator/(double a) const { 
		Matrix<_numRows, _numColumns> m;
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			m._elems[i] = _elems[i] / a;
		}
		return m;
	}
	inline const Matrix<_numRows, _numColumns>& operator/=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}

	// Matrix multiplication
	template <size_t numColumns>
	inline Matrix<_numRows, numColumns> operator*(const Matrix<_numColumns, numColumns>& q) const {
		Matrix<_numRows, numColumns> m;
		for (size_t i = 0; i < _numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				double temp = double(0);
				for (size_t k = 0; k < _numColumns; ++k) {
					temp += (*this)(i, k) * q(k, j);
				}
				m(i, j) = temp;
			}
		}
		return m;
	}

	inline const Matrix<_numRows, _numColumns>& operator*=(const Matrix<_numColumns, _numColumns>& q) { 
		return ((*this) = (*this) * q); 
	}

	// Matrix transpose
	inline Matrix<_numColumns, _numRows> operator~() const {
		Matrix<_numColumns, _numRows> m;
		for (size_t i = 0; i < _numColumns; ++i) {
			for (size_t j = 0; j < _numRows; ++j) {
				m(i, j) = (*this)(j, i);
			}
		}
		return m;
	}
};

/*template < >
class Matrix<1,1>: public Matrix
{
// Casting to double for 1x1 matrix
inline operator double() const {
return _elems[0];
}
};*/


// Scalar multiplication 
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(double a, const Matrix<_numRows, _numColumns>& q) { return q*a; }

// Matrix trace
template <size_t _size>
inline double tr(const Matrix<_size, _size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}

// Matrix 1-norm
template <size_t _numRows, size_t _numColumns>
inline double norm(const Matrix<_numRows, _numColumns>& q) {
	double norm1 = double(0);
	for (size_t j = 0; j < _numColumns; ++j) {
		double colabssum = double(0);
		for (size_t i = 0; i < _numRows; ++i) {
			colabssum += abs(q(i,j));
		}
		if (colabssum > norm1) {
			norm1 = colabssum;
		}
	}
	return norm1;
}


// Identity matrix
template <size_t _size>
inline Matrix<_size, _size> identity() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = (i == j ? double(1) : double(0));
		}
	}
	return m;
}

// Zero matrix
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> zeros() {
	Matrix<_numRows, _numColumns> m;
	m.reset();
	return m;
}

// Matrix determinant
template <size_t _size>
inline double det(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	double D = double(1);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = std::abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
			D = -D;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
			D = -D;
		}

		D *= m(row_p[k], col_p[k]);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		}  
	}

	return D;
}

// A%B solves PX = Q for X
template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator%(const Matrix<_size, _size>& p, const Matrix<_size, _numColumns>& q) {
	Matrix<_size, _size> m(p);
	Matrix<_size, _numColumns> inv(q);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		assert(maximum != double(0));

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);
		for (size_t j = 0; j < _numColumns; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// reshuffle result
	Matrix<_size, _numColumns> result;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			result(col_p[i], j) = inv(row_p[i], j);
		}
	}

	return result; 
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const Matrix<_size, _size>& q) {
	return ~(~q%~p);
}

// Matrix pseudo-inverse
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> pseudoInverse(const Matrix<_numRows, _numColumns>& q) { 
	Matrix<_numRows, _numColumns> m(q);

	size_t row_p[_numRows];
	size_t col_p[_numColumns];
	for (size_t i = 0; i < _numRows; ++i) {
		row_p[i] = i;
	}
	for (size_t i = 0; i < _numColumns; ++i) {
		col_p[i] = i;
	}

	// Gaussian elimination to determine rank
	size_t rank = (_numRows < _numColumns ? _numRows : _numColumns);
	for (size_t k = 0; k < rank; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _numRows; ++i) {
			for (size_t j = k; j < _numColumns; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// check if zero
		if (maximum <= DBL_EPSILON) {
			rank = k;
			break;
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
		}

		// eliminate column
		for (size_t i = k + 1; i < _numRows; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _numColumns; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		} 
	}

	if (rank == 0) { // zero or empty matrix
		return zeros<_numColumns, _numRows>();
	} else if (rank == _numRows && rank == _numColumns) { // full rank square matrix
		// convert so that compiler accepts number of rows and columns for inverse
		Matrix<_numRows, _numRows> q2;
		for (size_t i = 0; i < _numRows * _numRows; ++i) {
			q2[i] = q[i];
		}
		q2 = !q2;
		Matrix<_numColumns, _numRows> q3;
		for (size_t i = 0; i < _numRows * _numRows; ++i) {
			q3[i] = q2[i];
		}
		return q3;
	} else if (rank == _numRows) { // full row-rank matrix
		return ~q*!(q*~q);
	} else if (rank == _numColumns) { // full column-rank matrix
		return !(~q*q)*~q;
	} // else: not full rank, perform rank decomposition

	// bring m into reduced row echelon form

	// normalize rows such that 1's appear on the diagonal
	for (size_t k = 0; k < rank; ++k) {
		double quotient = m(row_p[k], col_p[k]);
		for (size_t j = 0; j < k; ++j) {
			m(row_p[k], col_p[j]) = double(0);
		}
		m(row_p[k], col_p[k]) = double(1);
		for (size_t j = k + 1; j < _numColumns; ++j) {
			m(row_p[k], col_p[j]) /= quotient;
		}
	}

	// subtract rows such that 0's appear above leading row elements
	for (size_t k = rank - 1; k != -1; --k) {
		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			m(row_p[i], col_p[k]) = double(0);
			for (size_t j = k + 1; j < _numColumns; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		} 
	}

	// copy m into smaller matrix C and swap columns
	Matrix<(_numRows < _numColumns ? _numRows : _numColumns), _numColumns> C;
	for (size_t k = 0; k < rank; ++k) {
		for (size_t j = 0; j < _numColumns; ++j) {
			C(k,j) = m(row_p[k], j);
		}
	}
	for (size_t k = rank; k < C.numRows(); ++k) {
		for (size_t j = 0; j < _numColumns; ++j) {
			C(k,j) = double(0);
		}
	}

	// B is copy of A with columns swapped and non-pivot columns left out
	Matrix<_numRows, (_numRows < _numColumns ? _numRows : _numColumns)> B;
	for (size_t k = 0; k < _numRows; ++k) {
		for (size_t j = 0; j < rank; ++j ) {
			B(k,j) = q(k, col_p[j]);
		}
		for (size_t j = rank; j < B.numColumns(); ++j) {
			B(k,j) = double(0);
		}
	}

	// Now, q = B*C and B and C are of full row/column rank
	return ~C*!((~B*B)*(C*~C))*~B;
}

// Matrix inverse
template <size_t _size>
inline Matrix<_size, _size> operator!(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	Matrix<_size, _size> inv = identity<_size>();

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
		}

		// eliminate column
		assert(m(row_p[k], col_p[k]) != double(0));
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
			for (size_t j = 0; j < k; ++j) {
				inv(row_p[i], row_p[j]) -= factor * inv(row_p[k], row_p[j]);
			}
			inv(row_p[i], row_p[k]) = -factor;
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);
		for (size_t j = 0; j < _size; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _size; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// reshuffle result
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(col_p[i], j) = inv(row_p[i], j);
		}
	}

	return m; 
}


template <size_t _size>
inline void jacobi(const Matrix<_size, _size>& q, Matrix<_size, _size>& V, Matrix<_size, _size>& D) {
	D = q;
	V = identity<_size>();

	while (true) {
		double maximum = 0; size_t max_row = 0; size_t max_col = 0;
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = i + 1; j < _size; ++j) {
				if (abs(D(i,j)) > maximum) {
					maximum = abs(D(i,j));
					max_row = i;
					max_col = j;
				}
			}
		}

		if (maximum <= DBL_EPSILON) {
			break;
		}

		double theta = (D(max_col, max_col) - D(max_row, max_row)) / (2 * D(max_row, max_col));
		double t = 1 / (abs(theta) + sqrt(theta*theta+1));
		if (theta < 0) t = -t;
		double c = 1 / sqrt(t*t+1); 
		double s = c*t;

		Matrix<_size, _size> R = identity<_size>();
		R(max_row,max_row) = c;
		R(max_col,max_col) = c;
		R(max_row,max_col) = s;
		R(max_col,max_row) = -s;

		// update D // 
		//std::cout << D << std::endl;
		//D = ~R * D * R;

		double temp1 = c*c*D(max_row, max_row) + s*s*D(max_col, max_col) - 2*c*s*D(max_row, max_col);
		double temp2 = s*s*D(max_row, max_row) + c*c*D(max_col, max_col) + 2*c*s*D(max_row, max_col);
		D(max_row, max_col) = 0;
		D(max_col, max_row) = 0;
		D(max_row, max_row) = temp1;
		D(max_col, max_col) = temp2;
		for (int j = 0; j < _size; ++j) {
			if ((j != max_row) && (j != max_col)) {
				temp1 = c * D(j, max_row) - s * D(j, max_col);
				temp2 = c * D(j, max_col) + s * D(j, max_row);
				D(j, max_row) = (D(max_row, j) = temp1);
				D(j, max_col) = (D(max_col, j) = temp2);
			}
		}
		//std::cout << D << std::endl << std::endl;


		V = V * R;
	} 
}


// Matrix exponentiation
#define _B0 1729728e1
#define _B1 864864e1
#define _B2 199584e1
#define _B3 2772e2
#define _B4 252e2
#define _B5 1512e0
#define _B6 56e0
#define _B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <size_t _size>
inline Matrix<_size, _size> exp(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> A(q);
	int s = (int) std::max(double(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E)); 

	A /= pow(2.0,s);
	Matrix<_size, _size> A2(A*A);
	Matrix<_size, _size> A4(A2*A2);
	Matrix<_size, _size> A6(A2*A4);
	Matrix<_size, _size> U( A*(A6*_B7 + A4*_B5 + A2*_B3 + identity<_size>()*_B1) );
	Matrix<_size, _size> V( A6*_B6 + A4*_B4 + A2*_B2 + identity<_size>()*_B0 ); 
	Matrix<_size, _size> R7 = !(V - U)*(V + U);

	for (int i = 0; i < s; ++i) {
		R7 *= R7;
	}
	return R7;
}

// Input stream
template <size_t _numRows, size_t _numColumns>
inline std::istream& operator>>(std::istream& is, Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			is >> q(i,j);
		}
	}
	return is;
}

// Output stream
template <size_t _numRows, size_t _numColumns>
inline std::ostream& operator<<(std::ostream& os, const Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			os << q(i,j) << "\t";
			//printf_s("%24.24g ", q(i,j));
		}
		//printf_s("\n");
		os << std::endl;
	}
	return os;
}

// Hilbert matrix
template <size_t _size>
inline Matrix<_size, _size> hilbert() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = double(1) / (double) (i + j + 1);
		}
	}
	return m;
}

// Solution to Continuous-time Algebraic Ricatti Equation ~AX + XA - XGX + Q = 0
template <size_t _size>
inline Matrix<_size, _size> care(const Matrix<_size, _size>& A, const Matrix<_size, _size>& G, const Matrix<_size, _size>& Q) {
	Matrix<2*_size, 2*_size> H, H2, W, SqrtH2;
	H.insert(0,0, A);      H.insert(0, _size, -G);
	H.insert(_size,0, -Q); H.insert(_size, _size, -~A);

	H2 = H*H;

	SqrtH2 = identity<2*_size>();
	for (size_t i = 0; i < 100; ++i) {
		SqrtH2 = 0.5*(SqrtH2 + !SqrtH2*H2);
	}

	W = H - SqrtH2;

	return W.subMatrix<_size, _size>(_size, 0) * !W.subMatrix<_size, _size>(0, 0);
}

// Eigen decomposition of a real symmetric matrix by converting to tridiagonal form followed by QL iteration
// From: Numerical recipes in C++ (3rd edition)
template <size_t _size>
inline void jacobi2(const Matrix<_size, _size>&	q, Matrix<_size, _size>& V, Matrix<_size, _size>& D) 
{
	// Initializations
	Matrix<_size, _size> z = q;
	Matrix<_size> d, e;
	d.reset();
	e.reset();

	int l,k,j,i,m,iter;
	double scale,hh,h,g,f,s,r,p,dd,c,b,absf,absg,pfg;

	//Convert to tridiagonal form
	for (i=_size-1;i>0;i--) 
	{
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z(i,k));
			if (scale == 0.0)
				e[i]=z(i,l);
			else {
				for (k=0;k<i;k++) {
					z(i,k) /= scale;
					h += z(i,k)*z(i,k);
				}
				f=z(i,l);
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z(i,l)=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					z(j,i)=z(i,j)/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z(j,k)*z(i,k);
					for (k=j+1;k<i;k++)
						g += z(k,j)*z(i,k);
					e[j]=g/h;
					f += e[j]*z(i,j);
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z(i,j);
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z(j,k) -= (f*e[k]+g*z(i,k));
				}
			}
		} else {
			e[i]=z(i,l);
		}
		d[i]=h;
	}
	d[0]=0.0;

	e[0]=0.0;
	for (i=0;i<_size;i++) {
		if (d[i] != 0.0) {
			for (j=0;j<i;j++) {
				g=0.0;
				for (k=0;k<i;k++)
					g += z(i,k)*z(k,j);
				for (k=0;k<i;k++)
					z(k,j) -= g*z(k,i);
			}
		}
		d[i]=z(i,i);
		z(i,i)=1.0;
		for (j=0;j<i;j++) z(j,i)=z(i,j)=0.0;
	}

	// QL iteration
	for (i=1;i<_size;i++) e[i-1]=e[i];
	e[_size-1]=0.0;
	for (l=0;l<_size;l++) {
		iter=0;
		do {
			for (m=l;m<_size-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= std::numeric_limits<double>::epsilon()*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) { 
					std::cerr << "Too many iterations in tqli" << std::endl;
					std::exit(-1);
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				absg=abs(g);
				r = ( (absg > 1.0) ? absg*sqrt(1.0+(1.0/absg)*(1.0/absg)) : sqrt(1.0+absg*absg));
				g=d[m]-d[l]+e[l]/(g+((g>=0.0) ? fabs(r):-fabs(r)));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					absf=abs(f);
					absg=abs(g);
					pfg = (absf > absg ? absf*sqrt(1.0+(absg/absf)*(absg/absf)) : (absg == 0.0 ? 0.0 : absg*sqrt(1.0+(absf/absg)*(absf/absg))));
					e[i+1]=(r=pfg); 
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<_size;k++) {
						f=z(k,i+1);
						z(k,i+1)=s*z(k,i)+c*f;
						z(k,i)=c*z(k,i)-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	
	// Sort eigenvalues and reorder eigenvectors
	/*
	for (int i=0;i<_size-1;i++) {
		double p=d[k=i];
		for (int j=i;j<_size;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (int j=0;j<_size;j++) {
				p=z(j,i);
				z(j,i)=z(j,k);
				z(j,k)=p;
			}
		}
	}
	*/

	//Copy over data
	V = z;
	D.reset();
	for(size_t i=0;i<_size;++i) {
		D(i,i) = d[i];
	}
}

template <size_t _size>
inline void forcePD(Matrix<_size, _size>& q) 
{
	Matrix<_size, _size> evec, eval;
	jacobi(q, evec, eval);
	for (size_t i = 0; i < _size; ++i) {
		if (eval(i,i) < 0) {
			eval(i,i) = DBL_EPSILON;
		}
	}
	q = evec * eval * ~evec;
}

template <size_t _numRows, size_t _numColumns>
inline void serializeMatrix(std::ofstream& ofs, const Matrix<_numRows, _numColumns>& q) 
{
	for (size_t i = 0; i < q.numRows(); ++i) {
		for (size_t j = 0; j < q.numColumns(); ++j) {
			ofs << std::setprecision(12) << " " << q(i,j);
		}
	}
}

template <size_t _numRows, size_t _numColumns>
inline void deserializeMatrix(std::ifstream& ifs, Matrix<_numRows, _numColumns>& q) 
{
	for (size_t i = 0; i < q.numRows(); ++i) {
		for (size_t j = 0; j < q.numColumns(); ++j) {
			ifs >> q(i,j);
		}
	}
}

template <size_t _numRows, size_t _numColumns>
inline void checkNaN(const Matrix<_numRows, _numColumns>& q) 
{
	for(size_t i = 0; i < q.numRows(); ++i) {
		for(size_t j = 0; j < q.numColumns(); ++j) {
			if (q(i,j) != q(i,j)) {
				std::cerr << "NaN values encountered! " << q << std::endl;
				std::exit(-1);
			}
		}
	}
}

// Given a positive-definite symmetric matrix A[0..n-1][0..n-1], construct and store its Cholesky decomposition, A = LLt.
template <size_t _size>
inline void cholesky(const Matrix<_size, _size>& A, Matrix<_size, _size>& L) 
{
	int n = (int)_size;
	L = A;
	double sum;
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			sum = L(i,j);
			for (int k = i-1; k >= 0; k--) sum -= L(i,k)*L(j,k);
			if (i == j) {
				if (sum <= 0) {
					std::cerr << "A, with rounding errors, is not positive-definite, cholesky failed" << std::endl;
					std::exit(-1);
				}
				L(i,i) = sqrt(sum);
			} 
			else {
				L(j,i) = sum/L(i,i);
			}
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			L(j,i) = L(i,j);
			L(i,j) = 0.0;
		}
	}
}

#endif