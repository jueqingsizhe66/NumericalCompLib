#ifndef __GAUSS_SEIDEL_SOR_ITERATION_H__
#define __GAUSS_SEIDEL_SOR_ITERATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Gauss-Seidel SOR iteration, basing on Jacobi iteration, is to solve large scale equation groups
// using iterative approach
// one important parameter is a factor, w
// Algorithm:
// (1) Given A, b(Ax = b) and x_0; k = 0, epsilon, N and a factor w
// (2) for 1<=i<=n  x_i(k+1) = (1 - w) * x_k + (w / a_ii) (b_i - SUM(1<=j<=i-1)(a_ij * x_j(k+1)) - SUM(i+1<=j<=n)(a_ij * x_j(k)))
// (3) if (||x(k+1) - x(k)||_infinite < epsilon) stop; return x(k+1);
// (4) if (k==N) stop; else ++k;

/* Includes */
#include "Jacobi_iteration.h"

/* Declaration */
namespace Felix { // namespace begins here

// this is a generator for factor w
template <typename A_T> 
typename A_T::ele_type Gauss_Seidel_SOR_factor_generator(const A_T &a) {
	return typename A_T::ele_type(2) / (typename A_T::ele_type(1) + pow((typename A_T::ele_type(1) + pow(a.get_spectral_radius(), 2)), typename A_T::ele_type(0.5)));
}
 
// Gauss-Seidel SOR iteration 
template <typename A_T, typename CNT_T, typename PRINT_T>
Jacobi_Gauss_Seidel_iteration_returnT<A_T> Gauss_Seidel_SOR_iteration(const A_T &a, const A_T &b, const typename A_T::ele_type &w, const typename A_T::ele_type &epsilon, const CNT_T &N, PRINT_T p) {
    CNT_T k(0);
	A_T x_k(a.get_row_size(), 1); 
	typename A_T::ele_type diff(epsilon + typename A_T::ele_type(1));

	while (k < N && diff > epsilon) {
		p << "k:  " << k << '\n';
		p << "x:\n" << x_k;
		(k == 0) ? (p << "diff:  ---\n\n") : (p << "diff:  " << diff << "\n\n");

		// use the same x sequence
		diff = Felix::Jacobi_Gauss_Seidel_step_iteration(a, b, x_k, x_k, w);
		++k;
	} 
	if (k == N)
		return Felix::Jacobi_Gauss_Seidel_iteration_returnT<A_T>(false, x_k);
	else
		return Felix::Jacobi_Gauss_Seidel_iteration_returnT<A_T>(true, x_k);
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[3][3] = { {4, 3, 0},
	                       {3, 4, -1},
	                       {0, -1, 4} };
	double b_arr[3][1] = { {24},
		                   {30},
						   {-24} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	double w = Felix::Gauss_Seidel_SOR_factor_generator(A);
	Felix::Jacobi_Gauss_Seidel_iteration_returnT<Felix::matrix<double> > res =
		Felix::Gauss_Seidel_SOR_iteration(A, b, w, 0.001, 30, Felix::Stream_Printer(std::cout));

	if (res.get_valid()) {
		std::cout << "w: " << w << std::endl;
		std::cout << res.get_x();
	}
	else {
		std::cerr << "Root not found.";
	}

	return 0;
}

*/

