#ifndef __STEFFENSEN_ACCELERATED_ITERATIVE_SEARCH_H__
#define __STEFFENSEN_ACCELERATED_ITERATIVE_SEARCH_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Steffensen method is an improved iterative method to solve non-linear equation
// Algorithm:
// Given f(x), x = f(x)
// (1) Given x_0, N, epsilon, and f(x).  k = 0;
// (2) y_k = f(x_k); z_k = f(y_k);
// (3) x_k+1 = x_k - (y_k - x_k) ^ 2 / (z_k - 2 * y_k + x_k);
// (4) if (|x_k+1 - x_k| < epsilon) stop;
// (5) if (k == N) stop; else ++k;
// Moreover, it can print all the computing procedures

/* Includes */
#include "iterative_search.h"

/* Declaration */
namespace Felix { // namespace begins here

template <typename FUNC_T, typename CNT_T, typename PRINT_T>
Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> Steffensen_accelerated_iterative_search(const FUNC_T &f, typename FUNC_T::x_type x_k, const CNT_T &N, 
										 const typename FUNC_T::x_type &epsilon, PRINT_T p) {
    CNT_T k(0);
	typename FUNC_T::x_type y_k;
	typename FUNC_T::x_type z_k;
	typename FUNC_T::x_type x_k_1(x_k);
	typename FUNC_T::x_type diff(epsilon + 1);

	p << "k\tx\ty\tz\n";
	while (k < N && (!(Felix::abs(diff) < epsilon))) {
	    y_k = f(x_k);
	    z_k = f(y_k);
		p << k << '\t' << x_k << '\t' << y_k << '\t' << z_k << '\n';
		x_k_1 = x_k - (((y_k - x_k) * (y_k - x_k)) / (z_k - 2 * y_k + x_k));
		diff = x_k_1 - x_k;
		x_k = x_k_1;
		++k;
	}
    if (k == N) return Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type>(false, x_k_1);
	else return Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type>(true, x_k_1);
}

} // namespace ends here


#endif

/* Test Case & User Manual

F should be q(x) where x = q(x).

class F {
public:
	typedef double x_type;
	typedef double y_type;

	double operator () (double x) const {
		return x * x * x  - 1;
	}
};

int main() {
	Felix::Non_linear_equations_search_returnT<double> res = 
		Felix::Steffensen_accelerated_iterative_search(Felix::make_function(F()), 1.5, 20, 0.001, Felix::Stream_Printer(std::cout));
	if (res.get_valid()) {
		std::cout << "Answer: " << res.get_x() << std::endl;
	}
	else {
		std::cerr << "Root not found." << std::endl;
	}

	return 0;
}

*/

