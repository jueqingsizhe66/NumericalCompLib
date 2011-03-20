#ifndef __NEWTON_BI_SEARCH_H__
#define __NEWTON_BI_SEARCH_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// This is a combination of Newton's method and bisection, which 
// makes use of the advantages of both.
// Algorithm:
// Given f(x), f(x) = 0
// First use bisection reduce the length of root interval to 1,
// then apply Newton method to converge quickly 

/* Includes */
#include "bisection.h"
#include "Newton_search.h"

/* Declaration */
namespace Felix { // namespace begins here
    
template <typename FUNC_T, typename FUNC_D_T, typename CNT_T, typename PRINT_T>
Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> Newton_bi_search(const FUNC_T &f, const FUNC_D_T &f_derivatives, typename FUNC_T::x_type a, typename FUNC_T::x_type b, 
								 const CNT_T &N, const typename FUNC_T::x_type &epsilon, PRINT_T p) {
    // first use bisection and set epsilon = 1
    Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> x_ = bisection(f, a, b, 1, p);     
 
	if (x_.get_valid())
        // then use Newton method
        return Newton_search(f, f_derivatives, x_.get_x(), N, epsilon, p);
	else
		return x_;
}


} // namespace ends here


#endif

/* Test Case & User Manual

F should be f(x) whose root you are searching for.
Fd = F'(x).

class F {
public:
	typedef double x_type;
	typedef double y_type;

	double operator () (double x) const {
		return x * x * x  - x - 1;
	}
};

class Fd {
public:
	typedef double x_type;
	typedef double y_type;

	double operator () (double x) const {
		return 3 * x * x - 1;
	}
};

int main() {
	Felix::Non_linear_equations_search_returnT<double> res = 
		Felix::Newton_bi_search(Felix::make_function(F()), Felix::make_function(Fd()), 0, 3.0, 20, 0.001, Felix::Stream_Printer(std::cout));
	if (res.get_valid()) {
		std::cout << "Answer: " << res.get_x() << std::endl;
	}
	else {
		std::cerr << "Root not found." << std::endl;
	}

	return 0;
}

*/

