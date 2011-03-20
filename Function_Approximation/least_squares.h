#ifndef __LEAST_SQUARES_H__
#define __LEAST_SQUARES_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Least squares can approximate a function using experimental data such as 
// (x_i, y_i) = (x_i, f(x_i))
// Algorithm:
// P_n(x) = a_0 + a_1*x + a_2*x*x + ... + a_n*x^n
// Given 0<=i<=m, x_i and y_i
// Ax = b
// A:                                                                          x:        =   b:
//  _                                                                     _    _     _       _                        _
// |  SUM(0<=i<=m)(1)      SUM(0<=i<=m)(x_i)    SUM(0<=i<=m)(x_i^2)   ...  |  |  a_0  |     |  SUM(0<=i<=m)(y_i)       |
// |  SUM(0<=i<=m)(x_i)    SUM(0<=i<=m)(x_i^2)  SUM(0<=i<=m)(x_i^3)   ...  |  |  a_1  |  =  |  SUM(0<=i<=m)(y_i*x_i)   |
// |  SUM(0<=i<=m)(x_i^2)  SUM(0<=i<=m)(x_i^3)  SUM(0<=i<=m)(x_i^4)   ...  |  |  a_2  |     |  SUM(0<=i<=m)(y_i*x_i^2) |
// |_ ...                  ...                  ...                   ... _|  |_ ... _|     |_ ...                    _|


/* Includes */
#include <utility>
using std::pair;

#include "../Basic_Components/matrix.h"
#include "../Basic_Components/polynomial.h"
#include "../Linear_Systems/Iterative_Methods/Gauss_Seidel_iteration.h"

/* Declaration */
namespace Felix { // namespace begins here

// return a polynomial
template <typename T>
polynomial<T> least_squares(const vector<pair<T, T> > &x_y_vec, size_t n) {
    matrix<T> A(n, n), b(n, 1);
	vector<T> exp_x_seq, exp_yx_seq;

	for (size_t k = 0; k <= T(2) * n - T(2); ++k) {
	    T sum(0), sum_yx(0); 
		for (size_t i = 0; i != x_y_vec.size(); ++i) {
		    T exp(1);
			for (size_t j = 0; j != k; ++j) {
			    exp *= x_y_vec[i].first;
			}
			sum += exp;
			if (k < n) {
			    sum_yx += exp * x_y_vec[i].second;
			}
		}
		exp_x_seq.push_back(sum);
		if (k < n) {
		    exp_yx_seq.push_back(sum_yx);
		}
	}

	// fill A and b
	for (size_t i = 0; i != n; ++i) {
		for (size_t j = 0; j != n; ++j) {
		    A[i][j] = exp_x_seq[i + j];
		}
	}
	for (size_t i = 0; i != n; ++i) {
	    b[i][0] = exp_yx_seq[i];
	}

	// solve Ax=b
	Jacobi_Gauss_Seidel_iteration_returnT<matrix<T> > res = Gauss_Seidel_iteration(A, b, matrix<T>::epsilon, RAND_MAX, Null_Printer());
	vector<T> coeff;
	for (size_t i = 0; i != res.get_x().get_row_size(); ++i)
		coeff.push_back(res.get_x()[i][0]);
	return polynomial<T>(coeff);
}


} // namespace ends here

#endif

/* Test Case & User Manual

int main() {
	double x[] = { 0.0, 0.2, 0.4, 0.6, 0.8 };
	double y[] = { 0.9, 1.9, 2.8, 3.3, 4.2 };
	std::vector<std::pair<double, double> > x_y;

	for (int i = 0; i != 5; ++i)
		x_y.push_back(std::make_pair(x[i], y[i]));

	Felix::polynomial<double> res = Felix::least_squares(x_y, 2);
	std::cout << res << std::endl;
	
	return 0;
}

*/

