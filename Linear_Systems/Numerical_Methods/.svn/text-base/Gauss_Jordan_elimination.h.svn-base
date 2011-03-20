#ifndef __GAUSS_JORDAN_ELIMINATION_H__
#define __GAUSS_JORDAN_ELIMINATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Gauss-Jordan elimination is a method to calculate the inverse matrix of A
// or solve equations
// Algorithm:
// (1) Given Ax = b, A, b.  det = 1
// (2) |a_rk| = max_k<=i<=n |aik|;
// (3) if (a_rk == 0) { det = 0; stop; }
// (4) if (r > k) { a_rk <=> a_rj; b_k <=> b_j; det *= -1; }
// (5) a_kj = a_kj / a_kk, k+1<=j<=n, b_k = b_k / a_kk;
//     det *= a_kk; a_kk = 1;
// (6) a_ij = a_ij - a_ik * a_kj, 1<=i<=n && i!=k, k+1<=j<=n; a_ik = 0;
// (7) b_i = b_i - a_ik * b_k, 1<=i<=n && i!=k

/* Includes */
#include "../../Basic_Components/Basic_Components.h"
#include "Gaussian_elimination.h"

/* Declaration */
namespace Felix { // namespace begins here

// This is an intrusive method to turn matrix a into E, and a's det
// Intrusive means the actual arguement will be altered in the function 
template <typename A_T, typename PRINT_T>
pair<bool, typename A_T::ele_type> intrusive_Gauss_Jordan_elimination(A_T &a, PRINT_T &p) {
	typename A_T::ele_type det(1);
	typename A_T::ele_type a_rk;
	size_t max_row;
	size_t k;
	for (k = 0; k < a.get_row_size(); ++k) {
		a_rk = a[k][k];
		max_row = k;
		p << a << '\n';
		for (size_t i = k + 1; i < a.get_row_size(); ++i) {
			if (Felix::abs(a_rk) < Felix::abs(a[i][k])) { // |a_rk| < |a_ik|
			    a_rk = a[i][k];
				max_row = i;
			}
		}
		if (Felix::abs(a_rk) < Felix::const_numbers_and_metrics<typename A_T::ele_type>::epsilon) {
			return make_pair(false, typename A_T::ele_type(0));
		}
		if (max_row > k) {
			a[k].swap(a[max_row]);
			det *= typename A_T::ele_type(-1);
		}

		for (size_t j = k + 1; j < a.get_col_size(); ++j)
			a[k][j] /= a[k][k];
		det *= a[k][k];
		a[k][k] = typename A_T::ele_type(1);
		for (size_t i = 0; i != a.get_row_size(); ++i) 
			if (i != k) {
			    for (size_t j = k + 1; j < a.get_col_size(); ++j)
					a[i][j] -= a[i][k] * a[k][j];
				a[i][k] = typename A_T::ele_type(0);
			}
	}
	p << a << '\n' << det << '\n';
	return make_pair(true, det);
}

// return type
// A -> bool & upper(A) & det(A)
template <typename A_T>
class Gauss_Jordan_elimination_partial_returnT { 
public:
	Gauss_Jordan_elimination_partial_returnT(bool v, const A_T &A, const typename A_T::ele_type &det) :
	  valid_(v), inverse_A_(A), det_(det) { }
	bool                   &get_valid() { return this->valid_; }
	A_T                    &get_inverse_A() { return this->inverse_A_; }
	typename A_T::ele_type &get_det() { return this->det_; }
private:
    bool                   valid_; // to check if the result is valid
    A_T                    inverse_A_;
	typename A_T::ele_type det_;
};

// A -> bool & upper(A) & det(A)
template <typename A_T, typename PRINT_T>
Gauss_Jordan_elimination_partial_returnT<A_T> Gauss_Jordan_elimination(A_T a, PRINT_T p) {
	pair<bool, typename A_T::ele_type> res = Felix::intrusive_Gauss_Jordan_elimination(a.horizontal_concat(A_T::get_E(a.get_row_size())), p);
	return Felix::Gauss_Jordan_elimination_partial_returnT<A_T>(res.first, a.get_local_matrix(0, a.get_row_size(), a.get_row_size(), a.get_col_size()), res.second);
}


// return type
// A * x[N] = b[N] -> bool & x[N] & inverse_A & det(A)
template <typename A_T, typename b_T = A_T>
class Gauss_Jordan_elimination_complete_returnT {
public:
	Gauss_Jordan_elimination_complete_returnT(bool v, const b_T &x, const A_T &A, const typename A_T::ele_type &det) :
	  valid_(v), x_(x), inverse_A_(A), det_(det) { }
	bool                   &get_valid() { return this->valid_; }
	b_T                    &get_x() { return this->x_; }
	A_T                    &get_inverse_A() { return this->inverse_A_; }
	typename A_T::ele_type &get_det() { return this->det_; }
private:
    bool                   valid_;
	b_T                    x_; // the root of the equations
    A_T                    inverse_A_;
	typename A_T::ele_type det_;
};

// typename A_T = b_T = matrix

/* Algorithm */
// After turn Ax = b into Ex = b', compute x: x = b'
// A * x[N] = b[N] -> bool & x[N] & inverse_A & det(A)
template <typename A_T, typename PRINT_T>
Gauss_Jordan_elimination_complete_returnT<A_T, vector<A_T> > Gauss_Jordan_elimination(A_T a, const vector<A_T> &b, PRINT_T p) {
	A_T b_chain(a.get_row_size(), 0);
	matrix_horizontal_concat(b_chain, b);
	a.horizontal_concat(b_chain);
	a.horizontal_concat(A_T::get_E(a.get_row_size()));
	// now matrix a is (a | b_0 | b_1 | ... | b_n-1 | E)
	pair<bool, typename A_T::ele_type> res = Felix::intrusive_Gauss_Jordan_elimination(a, p);
	// now matrix a is (E | x_0 | x_1 | ... | x_n-1 | inverse_a)
	if (res.first) {  
		vector<A_T> x(b.size(), A_T(a.get_row_size(), 1));
		for (size_t i = 0; i != b.size(); ++i) { // for each b[k] 
			for (size_t j = 0; j != a.get_row_size(); ++j)
				x[i][j][0] = a[j][a.get_row_size() + i];
			p << x[i] << '\n';
		}
		return Felix::Gauss_Jordan_elimination_complete_returnT<A_T, vector<A_T> >(true, x, a.get_local_matrix(0, a.get_row_size(), a.get_row_size() + b.size(), a.get_col_size()), res.second);
	}
	else // invalid
		return Felix::Gauss_Jordan_elimination_complete_returnT<A_T, vector<A_T> >(false, b, A_T(1, 1), typename A_T::ele_type(0));
}


// Ax = b -> bool & x & upper(A) & det(A)
template <typename A_T, typename PRINT_T>
Gauss_Jordan_elimination_complete_returnT<A_T, A_T> Gauss_Jordan_elimination(A_T a, A_T b, PRINT_T p) {
	vector<A_T> one_b; 
	one_b.push_back(b);
	Felix::Gauss_Jordan_elimination_complete_returnT<A_T, vector<A_T> > res = Felix::Gauss_Jordan_elimination(a, one_b, p);
	return Felix::Gauss_Jordan_elimination_complete_returnT<A_T, A_T>(res.get_valid(), (res.get_x())[0], res.get_inverse_A(), res.get_det());
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[3][3] = { {0, 1, 1},
	                       {1, 2, 3},
	                       {1, 4, 6} };
	double b_arr[3][1] = { {1},
		                   {0},
						   {2} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	Felix::Gauss_Jordan_elimination_partial_returnT<Felix::matrix<double> > res_partial = 
		Felix::Gauss_Jordan_elimination(A, Felix::Stream_Printer(std::cout));
	if (res_partial.get_valid()) {
		std::cout << "inverse(A):\n" << res_partial.get_inverse_A() << std::endl;
		std::cout << "det(A):\n" << res_partial.get_det() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	Felix::Gauss_Jordan_elimination_complete_returnT<Felix::matrix<double> > res_full = 
		Felix::Gauss_Jordan_elimination(A, b, Felix::Stream_Printer(std::cout));
	if (res_full.get_valid()) {
		std::cout << "x:\n" << res_full.get_x() << std::endl;
		std::cout << "inverse(A):\n" << res_full.get_inverse_A() << std::endl;
		std::cout << "det(A):\n" << res_full.get_det() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	std::vector<Felix::matrix<double> > b_vec;
	b_vec.push_back(b);
	Felix::Gauss_Jordan_elimination_complete_returnT<Felix::matrix<double>, std::vector<Felix::matrix<double> > > res_full_vec = 
		Felix::Gauss_Jordan_elimination(A, b_vec, Felix::Stream_Printer(std::cout));
	if (res_full_vec.get_valid()) {
		std::vector<Felix::matrix<double> > x = res_full_vec.get_x();
		std::cout << "x:\n" << x[0] << std::endl;
		std::cout << "inverse(A):\n" << res_full_vec.get_inverse_A() << std::endl;
		std::cout << "det(A):\n" << res_full_vec.get_det() << std::endl;
	}
	else {
	    std::cerr << "No valid result got.";
	}

	return 0;
}

*/
