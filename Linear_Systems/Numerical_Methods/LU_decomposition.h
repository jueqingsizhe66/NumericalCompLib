#ifndef __LU_DECOMPOSITION_H__
#define __LU_DECOMPOSITION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// LU decompostion is to turn A into A = LU; L is a lower triangular matrix and
// U is an upper triangular matrix
// By doing this, a more compacted form of Gaussian elimination is forged, then
// can be used to solve equation group.
// Algorithm:
// (1) Given A, b
// (2) for k=1 to n, u_kj = a_kj - SUM(1<=q<=k-1)(l_kq * u_qj), k<=j<=n;
// (3) l_ik = (a_ik - SUM(1<=q<=k-1)(l_iq * u_qk)) / u_kk, k+1<=i<=n;
// (4) solve Ly = b, for k=1 to n, y_k = b_k - SUM(1<=j<=k-1)(l_kj * y_j)
// (5) solve Ux = y, for k=n to 1, x_k = (y_k - SUM(k+1<=j<=n)(u_kj * x_j)) / u_kk

/* Includes */
#include "../../Basic_Components/matrix.h"
#include "../../Auxiliary_Components/Printer.h"

/* Declaration */
namespace Felix { // namespace begins here

// return type
// A -> bool & L & U
template <typename A_T>
class LU_decomposition_partial_returnT {
public:
	LU_decomposition_partial_returnT(bool v, const A_T &L, const A_T &U) : valid_(v), L_(L), U_(U) { }
	bool &get_valid() { return this->valid_; }
	A_T &get_L() { return this->L_; }
	A_T &get_U() { return this->U_; }
private:
    bool valid_;
	A_T L_;
	A_T U_;
};


// A -> bool & L & U
template <typename A_T, typename PRINT_T>
LU_decomposition_partial_returnT<A_T> LU_decomposition(const A_T &a, PRINT_T p) {
	A_T L(A_T::get_E(a.get_row_size())), U(a.get_row_size(), a.get_row_size());
	for (size_t k = 0; k != a.get_row_size(); ++k) {
		p << L << '\n' << U << '\n';
		for (size_t j = k; j < a.get_row_size(); ++j) {
			typename A_T::ele_type sum(0); 
		    for (size_t q = 0; q < k; ++q)
				sum += L[k][q] * U[q][j];
			U[k][j] = a[k][j] - sum;
		}
		for (size_t i = k + 1; i < a.get_row_size(); ++i) {
		    typename A_T::ele_type sum(0); 
			for (size_t q = 0; q < k; ++q)
				sum += L[i][q] * U[q][k];
			if (Felix::abs(U[k][k]) < Felix::const_numbers_and_metrics<typename A_T::ele_type>::epsilon) {
				p << L << '\n' << U << '\n';
				return Felix::LU_decomposition_partial_returnT<A_T>(false, L, U);
			}
			else
			    L[i][k] = (a[i][k] - sum) / U[k][k];
		}
	}
	p << L << '\n' << U << '\n';
	return Felix::LU_decomposition_partial_returnT<A_T>(true, L, U);
}


// return type
// A -> bool & x & L & U
template <typename A_T, typename b_T = A_T>
class LU_decomposition_complete_returnT {
public:
	LU_decomposition_complete_returnT(bool v, const b_T &x, const A_T &L, const A_T &U) : valid_(v), x_(x), L_(L), U_(U) { }
	bool &get_valid() { return this->valid_; }
	b_T &get_x() { return this->x_; }
	A_T &get_L() { return this->L_; }
	A_T &get_U() { return this->U_; }
private:
    bool valid_;
	b_T x_;
	A_T L_;
	A_T U_;
};

/* Algorithm */
// After turn Ax = b into L(Ux[N]) = b[N]
// solve L y[N] = b[N] -> y[N] = U x[N]
// solve x[N]
// A * x[N] = b[N] -> bool & x[N] & L & U
template <typename A_T, typename PRINT_T>
LU_decomposition_complete_returnT<A_T, vector<A_T> > LU_decomposition(const A_T &a, const vector<A_T> &b, PRINT_T p) {
	Felix::LU_decomposition_partial_returnT<A_T> res = Felix::LU_decomposition(a, p);
	if (res.get_valid()) { // non-singular matrix
		vector<A_T> x;
		A_T y(a.get_row_size(), 1), tmp_x(a.get_row_size(), 1);
		for (size_t i = 0; i != b.size(); ++i) { // for each b[N]
		    // solve L y = b
			for (size_t k = 0; k != a.get_row_size(); ++k) {
				typename A_T::ele_type sum(0);
				for (size_t j = 0; j < k; ++j)
					sum += res.get_L()[k][j] * y[j][0];
				y[k][0] = b[i][k][0] - sum;
			}
			p << y << '\n';
			// solve U x = y
			for (size_t rk = 0; rk != a.get_row_size(); ++rk) {
			    size_t k = a.get_row_size() - 1 - rk;
				typename A_T::ele_type sum(0);
				for (size_t j = k + 1; j < a.get_row_size(); ++j)
					sum += res.get_U()[k][j] * tmp_x[j][0];
				tmp_x[k][0] = (y[k][0] - sum) / res.get_U()[k][k];
			}
			p << tmp_x << '\n';
			// add tmp_x to x
			x.push_back(tmp_x);
		}
		return Felix::LU_decomposition_complete_returnT<A_T, vector<A_T> >(true, x, res.get_L(), res.get_U());
	}
	else 
		return Felix::LU_decomposition_complete_returnT<A_T, vector<A_T> >(false, b, res.get_L(), res.get_U());
}

// A * x = b -> bool & x & L & U
template <typename A_T, typename PRINT_T>
LU_decomposition_complete_returnT<A_T> LU_decomposition(const A_T &a, const A_T &b, PRINT_T p) {
	vector<A_T> one_b;
	one_b.push_back(b);
	Felix::LU_decomposition_complete_returnT<A_T, vector<A_T> > res = Felix::LU_decomposition(a, one_b, p);
	return Felix::LU_decomposition_complete_returnT<A_T, A_T>(res.get_valid(), (res.get_x())[0], res.get_L(), res.get_U());
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[3][3] = { {3, -1, 2},
	                       {1,  2, 3},
	                       {2, -2, -1} };
	double b_arr[3][1] = { {12},
		                   {11},
						   {2} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	Felix::LU_decomposition_partial_returnT<Felix::matrix<double> > res_partial =
		Felix::LU_decomposition(A, Felix::Stream_Printer(std::cout));

	if (res_partial.get_valid()) {
		std::cout << "L:\n" << res_partial.get_L() << std::endl;
		std::cout << "U:\n" << res_partial.get_U() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	Felix::LU_decomposition_complete_returnT<Felix::matrix<double> > res_full =
		Felix::LU_decomposition(A, b, Felix::Stream_Printer(std::cout));

	if (res_full.get_valid()) {
		std::cout << "x:\n" << res_full.get_x() << std::endl;
		std::cout << "L:\n" << res_full.get_L() << std::endl;
		std::cout << "U:\n" << res_full.get_U() << std::endl;
	}
	else {
		std::cerr << "No valid result got.";
	}

	std::vector<Felix::matrix<double> > b_vec;
	b_vec.push_back(b);
	Felix::LU_decomposition_complete_returnT<Felix::matrix<double>, std::vector<Felix::matrix<double> > > res_full_vec = 
		Felix::LU_decomposition(A, b_vec, Felix::Stream_Printer(std::cout));

	if (res_full_vec.get_valid()) {
		std::vector<Felix::matrix<double> > x = res_full_vec.get_x();
	    std::cout << "x:\n" << x[0] << std::endl;
		std::cout << "L:\n" << res_full.get_L() << std::endl;
		std::cout << "U:\n" << res_full.get_U() << std::endl;
	}
	else {
	    std::cerr << "No valid result got.";
	}

	return 0;
}

*/

