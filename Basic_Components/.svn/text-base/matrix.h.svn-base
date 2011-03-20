#ifndef __MATRIX_H__
#define __MATRIX_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Matrix is a very important tool in mathematical calculation
// This class can provide users with powerful matrix operations

/* Includes */
#include <iostream>
using std::istream;
using std::ostream;

#include <cmath>

#include <vector>
using std::vector;

#include "basic_widgets.h"
#include "../Linear_Systems/Numerical_Methods/Gaussian_elimination.h"
#include "../Linear_Systems/Numerical_Methods/Gauss_Jordan_elimination.h"

/* Declaration */
namespace Felix { // namespace Felix begins here.

// forward declaration
template <typename T>
class matrix;

template <typename T>
ostream &operator << (ostream &os, const matrix<T> &mtx);

template <typename T>
istream &operator >> (istream &is, matrix<T> &mtx);

template <typename T>
matrix<T> operator + (const matrix<T> &, const matrix<T> &);

template <typename T>
matrix<T> operator - (const matrix<T> &, const matrix<T> &);

template <typename T>
matrix<T> operator * (const matrix<T> &, const T &);

template <typename T>
matrix<T> operator * (const T &, const matrix<T> &);

template <typename T>
matrix<T> operator * (const matrix<T> &, const matrix<T> &);

template <typename T>
matrix<T> operator / (const matrix<T> &, const T &);

template <typename T>
bool operator == (const matrix<T> &lhs, const matrix<T> &rhs);

template <typename T>
bool operator != (const matrix<T> &lhs, const matrix<T> &rhs);


// Matrix Model
template <typename T> 
class matrix {
public:
	// define types
	typedef vector<vector<T> > con_type;
    typedef vector<T> row_type;
	typedef vector<T> col_type;
    typedef T         ele_type;
	// common interfaces
	// constructor
	matrix(size_t row, size_t col) : data_(row, vector<T>(col, 0)) { }

	template <size_t ROW, size_t COL>
	matrix(const T (&arr)[ROW][COL]) : data_(ROW) { for (size_t i = 0; i != ROW; ++i) this->data_[i].assign(arr[i], arr[i] + COL); }

	// copy constructor --- default
	// operator = --- default
	// destructor --- default

	size_t get_row_size() const { return this->data_.size(); }
	size_t get_col_size() const { return this->data_.size() == 0 ? 0 : this->data_[0].size(); }
	con_type &get_container() { return this->data_; }
	row_type &operator [] (size_t row) { if (row < this->data_.size()) return this->data_[row]; else throw computation_exception("OPERATOR [] OVERFLOWS!"); }
	const row_type &operator [] (size_t row) const { if (row < this->data_.size()) return this->data_[row]; else throw computation_exception("OPERATOR [] OVERFLOWS!"); }

	// matrix eigenvalues and eigenvectors ********************************************************
    // eigenvalue problem is so important for a matrix, so it deserves special treatment
	
	// return type
	class Eigenvalue_and_eigenvector {
	public:
		Eigenvalue_and_eigenvector(const T &v, const matrix<T> &m) : val_(v), vec_(m) { }
		T &get_eigenvalue() { return this->val_; }
		matrix<T> &get_eigenvector() { return this->vec_; }
	private:
		T val_;
		matrix<T> vec_;
	};

	typedef Eigenvalue_and_eigenvector eigenvalue_and_eigenvector_pair_type;
	typedef vector<eigenvalue_and_eigenvector_pair_type> eigenvalue_and_eigenvector_vector_type;

    // compute the eigenvalue and eigenvector with largest magnitude using power iteration
	eigenvalue_and_eigenvector_pair_type get_max_abs_eigenvalue_and_eigenvector_using_power_iteration() const;
	// compute the eigenvalue and eigenvector with smallest magnitude using inverse power iteration 
	eigenvalue_and_eigenvector_pair_type get_min_abs_eigenvalue_and_eigenvector_using_inverse_power_iteration() const;
	// compute an eigenvalue and eigenvector using an estimate of eigenvalue
	eigenvalue_and_eigenvector_pair_type get_an_eigenvalue_and_eigenvector_using_inverse_power_iteration(const T &) const;
	eigenvalue_and_eigenvector_pair_type get_an_eigenvalue_and_eigenvector_using_Rayleigh_quotient_iteration(T) const;


	// end of eigenvalue and eigenvector **********************************************************



	// properties of matrix
	// norms
	T get_Frobenius_norm() const;
	T get_1_norm() const;
	T get_2_norm() const; 
	T get_infinite_norm() const;
	// condition numbers
	T get_spectral_radius() const; 
	T get_1_condition_number() const;
	T get_1_condition_number(const matrix<T> &inverse_matrix) const;
	T get_2_condition_number() const; 
	T get_2_condition_number(const matrix<T> &inverse_matrix) const; 
	T get_infinite_condition_number() const; 
	T get_infinite_condition_number(const matrix<T> &inverse_matrix) const;
	// trace
	T get_trace() const;
	// special checks
	bool is_symmetrical() const;
	bool is_non_singular() const;

	// matrix operations
	matrix<T> &horizontal_concat(const matrix<T> &m);
	matrix<T> &vertical_concat(const matrix<T> &m);
	matrix<T> get_local_matrix(size_t row_begin, size_t row_end, size_t col_begin, size_t col_end);
	matrix<T> get_transpose_matrix() const;
	matrix<T> &transpose();

	// arithmetical operations
	friend matrix<T> operator + <T> (const matrix<T> &lhs, const matrix<T> &rhs);
	friend matrix<T> operator - <T> (const matrix<T> &lhs, const matrix<T> &rhs);
	friend matrix<T> operator * <T> (const matrix<T> &lhs, const matrix<T> &rhs);
	friend matrix<T> operator * <T> (const matrix<T> &lhs, const T &rhs);
	friend matrix<T> operator * <T> (const T &lhs, const matrix<T> &rhs);
	friend matrix<T> operator / <T> (const matrix<T> &lhs, const T &rhs);

	// relational operations
	friend bool operator == <T> (const matrix<T> &lhs, const matrix<T> &rhs);
	friend bool operator != <T> (const matrix<T> &lhs, const matrix<T> &rhs);

	// iostream operations
	friend ostream &operator << <T> (ostream &os, const matrix<T> &mtx);
	friend istream &operator >> <T> (istream &is, matrix<T> &mtx);

	// operations that have been implemented in other .h files, simply for convenience
	T get_determinant() const;
	matrix<T> get_inverse_matrix() const;
private:
	vector<vector<T> > data_;

	// static functions for special matrice
public:
    static matrix<T> get_E(size_t size);
	static matrix<T> get_Hilbert_ill_conditioned_matrix(size_t size);

	// precision unit
	static T epsilon;
};

/* Template Implementation */
// static member data
template <typename T>
T matrix<T>::epsilon = 0.0001;

// matrix eigenvalues and eigenvectors ********************************************************
// implementations of various eigenvalue problem solving algorithms

// compute the eigenvalue and eigenvector with largest magnitude using power iteration
// algorithm:
// (1) v_0 is random
// (2) for k = 1 to n
//     w = A * v_(k-1);
//     v_k = w / || w ||_infinite;
//     until no big change occurs
// (3) lamda = T(v)*A*v / T(v)*v is max eigenvalue and v is corresponding eigenvector    
template <typename T>
typename matrix<T>::eigenvalue_and_eigenvector_pair_type matrix<T>::get_max_abs_eigenvalue_and_eigenvector_using_power_iteration() const {
	matrix<T> v_k(this->get_row_size(), 1),
		      v_k_1(this->get_row_size(), 1),
		      w(this->get_row_size(), 1);
	// initiate v
	for (size_t i = 0; i != v_k_1.get_row_size(); ++i) {
		v_k_1[i][0] = typename matrix<T>::ele_type((Felix::positive_rand(0, 2) == 0 ? 1 : -1) * Felix::positive_rand(static_cast<T>(0), static_cast<T>(RAND_MAX)));
		v_k[i][0] = typename matrix<T>::ele_type(0);
	}

	while (Felix::abs((v_k_1 - v_k).get_infinite_norm()) > matrix<T>::epsilon) {
		v_k = v_k_1;
		w = (*this) * v_k_1;
		v_k_1 = w / w.get_infinite_norm();
	}
	return typename matrix<T>::eigenvalue_and_eigenvector_pair_type(
		   (v_k_1.get_transpose_matrix() * (*this) * v_k_1)[0][0] / (v_k_1.get_transpose_matrix() * v_k_1)[0][0], v_k_1);
}


// compute the eigenvalue and eigenvector with smallest magnitude using inverse power iteration
// simply call matrix<T>::get_an_eigenvalue_and_eigenvector_using_inverse_power_iteration with m = 0
template <typename T>
typename matrix<T>::eigenvalue_and_eigenvector_pair_type 
  matrix<T>::get_min_abs_eigenvalue_and_eigenvector_using_inverse_power_iteration() const {
    return this->get_an_eigenvalue_and_eigenvector_using_inverse_power_iteration(matrix<T>::ele_type(0));
}


// compute the eigenvalue and eigenvector using inverse power iteration
// Giving m to compute an eigenvalue near to m
// algorithm:
// (1) v_0 is random. Let A' = Inverse(A-mI), m should be as close to eigenvalue as possible
// (2) use Power Iteration to compute A' 's eigenvalue lamda' and eigenvector v'
// (3) lamda = 1/lamda' + m is max eigenvalue and v=v' is corresponding eigenvector    
template <typename T>
typename matrix<T>::eigenvalue_and_eigenvector_pair_type 
  matrix<T>::get_an_eigenvalue_and_eigenvector_using_inverse_power_iteration(const T &m) const {
	matrix<T> A_p = ((*this) - m * matrix<T>::get_E(this->get_row_size())).get_inverse_matrix();
	typename matrix<T>::eigenvalue_and_eigenvector_pair_type res = A_p.get_max_abs_eigenvalue_and_eigenvector_using_power_iteration();
	return matrix<T>::eigenvalue_and_eigenvector_pair_type(matrix<T>::ele_type(1) / res.get_eigenvalue() + m, res.get_eigenvector());
}


// compute the eigenvalue and eigenvector using Rayleigh quotient iteration
// Giving m to compute an eigenvalue
// It is a natural combination of power iteration and inverse power iteration
// It uses power iteration to estimate a m for inverse power iteration
// algorithm:
// (1) v_0 is random, m_0 = 0. Let A' = Inverse(A-m_k*I)
// (2) w = A' * v
// (3) v = w / ||w||
// (4) m_k+1 = T(v)*A*v / T(v)*v
// (5) goto (1) until no big change occurs to v
// (6) m_k+1 is max eigenvalue and v=v' is corresponding eigenvector    
template <typename T>
typename matrix<T>::eigenvalue_and_eigenvector_pair_type 
  matrix<T>::get_an_eigenvalue_and_eigenvector_using_Rayleigh_quotient_iteration(T m) const {
	matrix<T> v(this->get_row_size(), 1),
		      w(this->get_row_size(), 1),
	          A_p(this->get_row_size(), this->get_col_size());
	// initiate v
	for (size_t i = 0; i != v.get_row_size(); ++i) {
		v[i][0] = matrix<T>::ele_type((Felix::positive_rand(0, 2) == 0 ? 1 : -1) * Felix::positive_rand(static_cast<T>(0), static_cast<T>(RAND_MAX)));
	}

	// check if Av = mv
	while (Felix::abs(((*this) * v)[0][0] - (m * v)[0][0]) > matrix<T>::epsilon) {
		A_p = ((*this) - m * matrix<T>::get_E(this->get_row_size())).get_inverse_matrix();
		w = A_p * v;
		v = w / w.get_infinite_norm();
		m = (v.get_transpose_matrix() * (*this) * v)[0][0] / (v.get_transpose_matrix() * v)[0][0];
	}
	return matrix<T>::eigenvalue_and_eigenvector_pair_type(m, v);
}


// end of eigenvalue and eigenvector **********************************************************


// properties of matrix
// || A ||(F) = (SUM_SUM(a_ij ^ 2)) ^ 0.5
template <typename T>
T matrix<T>::get_Frobenius_norm() const {
	typename matrix<T>::ele_type sum(0);
	for (size_t i = 0; i != this->get_row_size(); ++i)
		for (size_t j = 0; j != this->get_col_size(); ++j)
			sum += this->data_[i][j] * this->data_[i][j];
	return pow(sum, 0.5);
}


// || A ||(1) = max[1<=i<=n](SUM_1<=j<=n a[j][i])
template <typename T>
T matrix<T>::get_1_norm() const {
	typename matrix<T>::ele_type sum(0),
		                max_sum(0);
	for (size_t i = 0; i != this->get_col_size(); ++i) {
		sum = matrix<T>::ele_type(0);
		for (size_t j = 0; j != this->get_row_size(); ++j) {
			sum += Felix::abs(this->data_[j][i]);
		}
		if (sum > max_sum)
			max_sum = sum;
	}
	return max_sum;
}


// || A ||(2) = lamda_max(Transpose_A * A) ^ 0.5
template <typename T>
T matrix<T>::get_2_norm() const {
	return pow(((this->get_transpose_matrix() * (*this)).get_max_abs_eigenvalue_and_eigenvector_using_power_iteration()).get_eigenvalue(), 0.5);
}


// || A ||(infinite) = max[1<=i<=n](SUM_1<=j<=n a[i][j])
template <typename T>
T matrix<T>::get_infinite_norm() const {
	typename matrix<T>::ele_type sum(0),
		                max_sum(0);
	for (size_t i = 0; i != this->get_row_size(); ++i) {
		sum = typename matrix<T>::ele_type(0);
		for (size_t j = 0; j != this->get_col_size(); ++j) {
			sum += Felix::abs(this->data_[i][j]);
		}
		if (sum > max_sum)
			max_sum = sum;
	}
	return max_sum;
}


// Spectral_radius(A) = |lamda|_max(Transpose_A * A)
template <typename T>
T matrix<T>::get_spectral_radius() const {
	return (this->get_max_abs_eigenvalue_and_eigenvector_using_power_iteration()).get_eigenvalue();
}


// cond(A) = || A || * || inverse(A) ||
// 1_condition_number
template <typename T>
T matrix<T>::get_1_condition_number() const {
	Felix::Gauss_Jordan_elimination_partial_returnT<matrix<T> > res = Felix::Gauss_Jordan_elimination(*this, Felix::Null_Printer());
	if (res.get_valid())
	    return this->get_1_norm() * res.get_inverse_A().get_1_norm();
	else
		throw computation_exception("CONDITION NUMBER INVALID FOR SINGULAR MATRIX!");
}

// if the user can provide inverse_matrix
template <typename T>
T matrix<T>::get_1_condition_number(const matrix<T> &inverse_matrix) const {
	return this->get_1_norm() * inverse_matrix.get_1_norm();
}

// 2_condition_number
template <typename T>
T matrix<T>::get_2_condition_number() const {
	Felix::Gauss_Jordan_elimination_partial_returnT<matrix<T> > res = Felix::Gauss_Jordan_elimination(*this, Felix::Null_Printer());
	if (res.get_valid())
	    return this->get_2_norm() * res.get_inverse_A().get_2_norm();
	else
		throw computation_exception("CONDITION NUMBER INVALID FOR SINGULAR MATRIX!");
}

// if the user can provide inverse_matrix
template <typename T>
T matrix<T>::get_2_condition_number(const matrix<T> &inverse_matrix) const {
	return this->get_2_norm() * inverse_matrix.get_2_norm();
}

// infinite_condition_number
template <typename T>
T matrix<T>::get_infinite_condition_number() const {
	Felix::Gauss_Jordan_elimination_partial_returnT<matrix<T> > res = Felix::Gauss_Jordan_elimination(*this, Felix::Null_Printer());
	if (res.get_valid())
	    return this->get_infinite_norm() * res.get_inverse_A().get_infinite_norm();
	else
		throw computation_exception("CONDITION NUMBER INVALID FOR SINGULAR MATRIX!");
}

// if the user can provide inverse_matrix
template <typename T>
T matrix<T>::get_infinite_condition_number(const matrix<T> &inverse_matrix) const {
	return this->get_infinite_norm() * inverse_matrix.get_infinite_norm();
}

// trace A: sum(a_ii) 1<=i<=n
template <typename T>
T matrix<T>::get_trace() const {
	typename matrix<T>::ele_type tr(0);
	for (size_t i = 0; i != this->get_row_size(); ++i)
		tr += (*this)[i][i];
	return tr;
}

// check if (*this) is a symmetrical matrix
template <typename T>
bool matrix<T>::is_symmetrical() const {
    if (this->get_row_size() != this->get_col_size())
		return false;
	for (size_t i = 0; i != this->get_row_size(); ++i)
		for (size_t j = 0; j != this->get_col_size(); ++j)
			if (Felix::abs((*this)[i][j] - (*this)[j][i]) > matrix<T>::epsilon)
				return false;
	return true;
}

// check if |A| != 0
template <typename T>
bool matrix<T>::is_non_singular() const {
    if (this->get_row_size() != this->get_col_size())
		return false;
	Felix::Gaussian_elimination_partial_returnT<matrix<T> > res = Felix::Gaussian_elimination((*this), Felix::Null_Printer());
	if (Felix::abs(res.get_det()) > matrix<T>::epsilon) // != 0
		return true;
	else
		return false;
}


// concatenate two matrice together horizontally
template <typename T>
matrix<T> &matrix<T>::horizontal_concat(const matrix<T> &m) {
    if (this->get_row_size() != m.get_row_size())
		throw computation_exception("MATRIX CONCATENATION ROW SIZE NOT MATCHED!");
	for (size_t i = 0; i != this->data_.size(); ++i) {
	    this->data_[i].insert((this->data_[i]).end(), (m.data_[i]).begin(), (m.data_[i]).end());
	}
	return *this;
}

// concatenate two matrice together vertically
template <typename T>
matrix<T> &matrix<T>::vertical_concat(const matrix<T> &m) {
    if (this->get_col_size() != m.get_col_size())
        throw computation_exception("MATRIX CONCATENATION COLUMN SIZE NOT MATCHED!");
	this->data_.insert(this->data_.end(), m.data_.begin(), m.data_.end());
	return *this;
}

// choose a part of matrix
template <typename T>
matrix<T> matrix<T>::get_local_matrix(size_t row_begin, size_t row_end, size_t col_begin, size_t col_end) {
	// part too big
	if (row_end > this->get_row_size() || col_end > this->get_col_size())
		return *this;
    matrix<T> res(0, 0);
	for (size_t i = row_begin; i < row_end; ++i)
		res.data_.push_back(vector<T>(this->data_[i].begin() + col_begin, this->data_[i].begin() + col_end));
	return res;
}

// matrix transposition
// non-intrusive
template <typename T>
matrix<T> matrix<T>::get_transpose_matrix() const {
	matrix<T> res(this->get_col_size(), this->get_row_size());
	for (size_t i = 0; i != this->get_row_size(); ++i)
		for (size_t j = 0; j != this->get_col_size(); ++j)
			res[j][i] = (*this)[i][j];
	return res;
}

// intrusive
template <typename T>
matrix<T> &matrix<T>::transpose() {
	return ((*this) = this->get_transpose_matrix());
}


// Arithmetical Operations
// addition
template <typename T>
matrix<T> operator + (const matrix<T> &lhs, const matrix<T> &rhs) {
    if ((lhs.get_row_size() != rhs.get_row_size()) || (lhs.get_col_size() != rhs.get_col_size()))
		throw computation_exception("MATRIX ADDITION SIZE NOT MATCHED!");
	matrix<T> sum(lhs.get_row_size(), lhs.get_col_size());
	for (size_t i = 0; i != lhs.get_row_size(); ++i)
		for (size_t j = 0; j != lhs.get_col_size(); ++j)
			sum[i][j] = lhs[i][j] + rhs[i][j];
	return sum;
}

// subtraction
template <typename T>
matrix<T> operator - (const matrix<T> &lhs, const matrix<T> &rhs) {
	if ((lhs.get_row_size() != rhs.get_row_size()) || (lhs.get_col_size() != rhs.get_col_size()))
		throw computation_exception("MATRIX SUBTRACTION SIZE NOT MATCHED!");
	matrix<T> diff(lhs.get_row_size(), lhs.get_col_size());
	for (size_t i = 0; i != lhs.get_row_size(); ++i)
		for (size_t j = 0; j != lhs.get_col_size(); ++j)
			diff[i][j] = lhs[i][j] - rhs[i][j];
	return diff;
}

// multiplication
template <typename T>
matrix<T> operator * (const matrix<T> &lhs, const matrix<T> &rhs) {
	if (lhs.get_col_size() != rhs.get_row_size())
		throw computation_exception("MATRIX MULTIPLICATION SIZE NOT MATCHED!");
	matrix<T> prod(lhs.get_row_size(), rhs.get_col_size());
	for (size_t i = 0; i != lhs.get_row_size(); ++i) 
		for (size_t j = 0; j != rhs.get_col_size(); ++j) 
			for (size_t k = 0; k != lhs.get_col_size(); ++k)
				prod[i][j] += lhs[i][k] * rhs[k][j];
	return prod;
}

template <typename T>
matrix<T> operator * (const matrix<T> &lhs, const T &rhs) {
    matrix<T> prod(lhs.get_row_size(), lhs.get_col_size());
	for (size_t i = 0; i != lhs.get_row_size(); ++i) 
		for (size_t j = 0; j != lhs.get_col_size(); ++j) 
			prod[i][j] = lhs[i][j] * rhs;
	return prod;
}

template <typename T>
matrix<T> operator * (const T &lhs, const matrix<T> &rhs) {
    matrix<T> prod(rhs.get_row_size(), rhs.get_col_size());
	for (size_t i = 0; i != rhs.get_row_size(); ++i) 
		for (size_t j = 0; j != rhs.get_col_size(); ++j) 
			prod[i][j] = lhs * rhs[i][j];
	return prod;
}

// division
template <typename T>
matrix<T> operator / (const matrix<T> &lhs, const T &rhs) {
    matrix<T> quot(lhs.get_row_size(), lhs.get_col_size());
	for (size_t i = 0; i != lhs.get_row_size(); ++i) 
		for (size_t j = 0; j != lhs.get_col_size(); ++j) 
			quot[i][j] = lhs[i][j] / rhs;
	return quot;
}


// Relational Operations
template <typename T>
bool operator == (const matrix<T> &lhs, const matrix<T> &rhs) {
	if ((lhs.get_row_size() != rhs.get_row_size()) || (lhs.get_col_size() != rhs.get_col_size()))
		return false;
	for (size_t i = 0; i != lhs.get_row_size(); ++i)
		for (size_t j = 0; j != lhs.get_col_size(); ++j)
			if (Felix::abs(lhs[i][j] - rhs[i][j]) > matrix<T>::epsilon)
				return false;
	return true;
}

template <typename T>
bool operator != (const matrix<T> &lhs, const matrix<T> &rhs) {
    return !(lhs == rhs);
}


// IO Operations
template <typename T>
ostream &operator << (ostream &os, const matrix<T> &mtx) {
	for (size_t i = 0; i != mtx.get_row_size(); ++i) {
		for (size_t j = 0; j != mtx.get_col_size(); ++j)
			os << mtx[i][j] << "\t";
		os << std::endl;
	}
	return os;
}

template <typename T>
istream &operator >> (istream &is, matrix<T> &mtx) {
	for (size_t i = 0; i != mtx.get_row_size(); ++i)
        for (size_t j = 0; j != mtx.get_col_size(); ++j)
			is >> mtx[i][j];
	return is;
}

template <typename T>
T matrix<T>::get_determinant() const {
	Felix::Gaussian_elimination_partial_returnT<matrix<T> > res = Felix::Gaussian_elimination(*this, Felix::Null_Printer());
	return res.get_det();
}

template <typename T>
matrix<T> matrix<T>::get_inverse_matrix() const {
	Felix::Gauss_Jordan_elimination_partial_returnT<matrix<T> > res = Felix::Gauss_Jordan_elimination(*this, Felix::Null_Printer());
	return res.get_inverse_A();
}


/* Static Functions for Special Matrices */
// unit matrix E
template <typename T>
matrix<T> matrix<T>::get_E(size_t size) {
    matrix<T> unit_matrix(size, size);
	for (size_t i = 0; i != size; ++i)
		unit_matrix[i][i] = T(1);
	return unit_matrix;
}

// Hilbert matrix
template <typename T>
matrix<T> matrix<T>::get_Hilbert_ill_conditioned_matrix(size_t size) {
    matrix<T> Hilbert_matrix(size, size);
	
	for (size_t i = 0; i != size; ++i)
		for (size_t j = 0; j != size; ++j)
			Hilbert_matrix[i][j] = matrix<T>::ele_type(matrix<T>::ele_type(1) / (i + j + matrix<T>::ele_type(1)));

	return Hilbert_matrix;
}


/* Convenient Matrix Operations */
// concatenate a chain of matrice together
template <typename A_T>
void matrix_horizontal_concat(A_T &ori, const vector<A_T> &vec) {
    for (size_t i = 0; i != vec.size(); ++i)
		ori.horizontal_concat(vec[i]);
}

} // namespace ends here

#endif
