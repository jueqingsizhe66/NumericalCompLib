#ifndef __GAUSS_SEIDEL_ITERATION_H__
#define __GAUSS_SEIDEL_ITERATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Gauss-Seidel iteration is a special case for Gauss_Seidel_SOR_iteration,
// where factor w equals 1


/* Includes */
#include "Gauss_Seidel_SOR_iteration.h"

/* Declaration */
namespace Felix { // namespace begins here


// Gauss-Seidel iteration 
template <typename A_T, typename CNT_T, typename PRINT_T>
Jacobi_Gauss_Seidel_iteration_returnT<A_T> Gauss_Seidel_iteration(const A_T &a, const A_T &b, const typename A_T::ele_type &epsilon, const CNT_T &N, PRINT_T p) {
	return Felix::Gauss_Seidel_SOR_iteration(a, b, typename A_T::ele_type(1), epsilon, N, p);
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[3][3] = { {11, -3,  -2},
	                       {-1,  5,  -3},
	                       {-2, -12, 19} };
	double b_arr[3][1] = { {3},
		                   {6},
						   {-7} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	Felix::Jacobi_Gauss_Seidel_iteration_returnT<Felix::matrix<double> > res =
		Felix::Gauss_Seidel_iteration(A, b, 0.001, 50, Felix::Stream_Printer(std::cout));

	if (res.get_valid()) {
		std::cout << res.get_x();
	}
	else {
		std::cerr << "Root not found.";
	}

	return 0;
}

*/

