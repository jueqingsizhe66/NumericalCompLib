#ifndef __GAUSSIAN_ELIMINATION_H__
#define __GAUSSIAN_ELIMINATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Gauss elimination is the simplest method to solve matrix problems
// or solve equations
// Algorithm:
// (1) Given Ax = b, A, b.  det = 1
// (2) |a_rk| = max_k<=i<=n |aik|;
// (3) if (a_rk == 0) { det = 0; stop; }
// (4) if (r > k) { a_rk <=> a_rj; b_k <=> b_j; det *= -1; }
// (5) m_ik = a_ik / a_kk; a_ij = a_ij - m_ik * a_kj; b_i = b_i - m_ik * b_k;
// (6) det *= a_kk;
// (7) if (a_nn == 0) { det = 0; stop; }
//     else x_k = b_k - sum_j=k+1->n(a_kj * x_j) / a_kk;
// (8) det *= a_nn * det;

/* Includes */
#include <algorithm>
#include <vector>
using std::vector;
#include <utility>
using std::pair;
using std::make_pair;

#include "../../Basic_Components/Basic_Components.h"
#include "../../Auxiliary_Components/Printer.h"

/* Declaration */
namespace Felix { // namespace begins here

// This is an intrusive method to compute upper triangular matrix of a, and its det
// Intrusive means the actual arguement will be altered in the function 
template <typename A_T, typename PRINT_T>
pair<bool, typename A_T::ele_type> intrusive_Gaussian_elimination(A_T &a, PRINT_T &p) {
    typename A_T::ele_type det(1);
	typename A_T::ele_type a_rk;
	size_t max_row;
	size_t k;
	for (k = 0; k < a.get_row_size() - 1; ++k) {
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

		for (size_t i = k + 1; i < a.get_row_size(); ++i) {
		    typename A_T::ele_type m_ik = a[i][k] / a[k][k];
			a[i][k] = typename A_T::ele_type(0);
			for (size_t j = k + 1; j < a.get_col_size(); ++j) 
				a[i][j] -= m_ik * a[k][j];
		} 
		det *= a[k][k];
	}
	if (Felix::abs(a[k][k]) < Felix::const_numbers_and_metrics<typename A_T::ele_type>::epsilon) {
		p << a << '\n' << typename A_T::ele_type(0) << '\n';
		return make_pair(false, typename A_T::ele_type(0));
	}
	else {
		det *= a[k][k];
		p << a << '\n' << det << '\n';
		return make_pair(true, det);
	}
}


// return type
template <typename A_T>
class Gaussian_elimination_partial_returnT { 
public:
	Gaussian_elimination_partial_returnT(bool v, const A_T &A, const typename A_T::ele_type &det) :
	  valid_(v), upper_A_(A), det_(det) { }
	bool                   &get_valid() { return this->valid_; }
	A_T                    &get_upper_A() { return this->upper_A_; }
	typename A_T::ele_type &get_det() { return this->det_; }
private:
    bool                   valid_; // to check if the result is valid
    A_T                    upper_A_;
	typename A_T::ele_type det_;
};

// A -> bool & upper(A) & det(A)
template <typename A_T, typename PRINT_T>
Gaussian_elimination_partial_returnT<A_T> Gaussian_elimination(A_T a, PRINT_T p) {
	pair<bool, typename A_T::ele_type> ans = Felix::intrusive_Gaussian_elimination(a, p);
	return Felix::Gaussian_elimination_partial_returnT<A_T>(ans.first, a, ans.second);
}


// return type
template <typename A_T, typename b_T = A_T>
class Gaussian_elimination_complete_returnT {
public:
	Gaussian_elimination_complete_returnT(bool v, const b_T &x, const A_T &A, const typename A_T::ele_type &det) :
	  valid_(v), x_(x), upper_A_(A), det_(det) { }
	bool                   &get_valid() { return this->valid_; }
	b_T                    &get_x() { return this->x_; }
	A_T                    &get_upper_A() { return this->upper_A_; }
	typename A_T::ele_type &get_det() { return this->det_; }
private:
    bool                   valid_;
	b_T                    x_; // the root of the equations
    A_T                    upper_A_;
	typename A_T::ele_type det_;
};

// typename A_T = b_T = matrix

/* Algorithm */
// After turn Ax = b into upper(A)x = b, compute x
// x_k = (b_k - SUM(j=k+1, n)(a_kj*x_j)) / a_kk;
// A * x[N] = b[N] -> bool & x[N] & upper(A) & det(A)
template <typename A_T, typename PRINT_T>
Gaussian_elimination_complete_returnT<A_T, vector<A_T> > Gaussian_elimination(A_T a, const vector<A_T> &b, PRINT_T p) {
	A_T b_chain(a.get_row_size(), 0);
	matrix_horizontal_concat(b_chain, b);
	a.horizontal_concat(b_chain);
	// now matrix a is (a | b_0 | b_1 | ... | b_n-1)
	pair<bool, typename A_T::ele_type> res = Felix::intrusive_Gaussian_elimination(a, p);
	if (res.first) {  
		vector<A_T> x(b.size(), A_T(a.get_row_size(), 1));
		for (size_t i = 0; i != b.size(); ++i) { // for each b[k] 
			for (size_t r_k = 0; r_k != a.get_row_size(); ++r_k) {
				size_t k = a.get_row_size() - 1 - r_k;
				typename A_T::ele_type sum(0);
				for (size_t j = k + 1; j < a.get_row_size(); ++j)
					sum += a[k][j] * x[i][j][0];
				x[i][k][0] = (a[k][a.get_row_size() + i] - sum) / a[k][k];
			}
			p << x[i] << '\n';
		}
		return Felix::Gaussian_elimination_complete_returnT<A_T, vector<A_T> >(true, x, a.get_local_matrix(0, a.get_row_size(), 0, a.get_row_size()), res.second);
	}
	else // invalid
		return Felix::Gaussian_elimination_complete_returnT<A_T, vector<A_T> >(false, b, A_T(1, 1), typename A_T::ele_type(0));
}

// Ax = b -> bool & x & upper(A) & det(A)
template <typename A_T, typename PRINT_T>
Gaussian_elimination_complete_returnT<A_T, A_T> Gaussian_elimination(A_T a, A_T b, PRINT_T p) {
	vector<A_T> one_b; 
	one_b.push_back(b);
	Felix::Gaussian_elimination_complete_returnT<A_T, vector<A_T> > res = Felix::Gaussian_elimination(a, one_b, p);
	return Felix::Gaussian_elimination_complete_returnT<A_T, A_T>(res.get_valid(), (res.get_x())[0], res.get_upper_A(), res.get_det());
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[2][2] = { {0.002, 87.13},
	                       {4.453, -7.26} };
	double b_arr[2][1] = { {87.15},
	                       {37.27} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	Felix::Gaussian_elimination_partial_returnT<Felix::matrix<double> > res_partial = 
		Felix::Gaussian_elimination(A, Felix::Stream_Printer(std::cout));
	if (res_partial.get_valid()) {
		std::cout << "upper(A):\n" << res_partial.get_upper_A() << std::endl;
		std::cout << "det(A):\n" << res_partial.get_det() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	Felix::Gaussian_elimination_complete_returnT<Felix::matrix<double> > res_full = 
		Felix::Gaussian_elimination(A, b, Felix::Stream_Printer(std::cout));
	if (res_full.get_valid()) {
		std::cout << "x:\n" << res_full.get_x() << std::endl;
		std::cout << "upper(A):\n" << res_full.get_upper_A() << std::endl;
		std::cout << "det(A):\n" << res_full.get_det() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	std::vector<Felix::matrix<double> > b_vec;
	b_vec.push_back(b);
	Felix::Gaussian_elimination_complete_returnT<Felix::matrix<double>, std::vector<Felix::matrix<double> > > res_full_vec = 
		Felix::Gaussian_elimination(A, b_vec, Felix::Stream_Printer(std::cout));
	if (res_full_vec.get_valid()) {
		std::vector<Felix::matrix<double> > x = res_full_vec.get_x();
		std::cout << "x:\n" << x[0] << std::endl;
		std::cout << "upper(A):\n" << res_full_vec.get_upper_A() << std::endl;
		std::cout << "det(A):\n" << res_full_vec.get_det() << std::endl;
	}
	else {
	    std::cerr << "No valid result got.";
	}

	return 0;
}

*/

