#ifndef __BISECTION_H__
#define __BISECTION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Bisection is an effective method in solving non-linear equation
// Algorithm:
// Given f(x), f(x) = 0
// (1) Given [a, b], f(x), and epsilon. Ensure f(a) * f(b) < 0
// (2) if ((b - a) / 2 < epsilon) stop
// (3) x = (a + b) / 2, if (f(a) * f(x) < 0) b = x; else a = x;
// (4) goto(2);
// Moreover, it can print all the computing procedures

/* Includes */
#include "iterative_search.h"

/* Declaration */
namespace Felix { // namespace begins here


template <typename FUNC_T, typename PRINT_T>
Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> bisection(const FUNC_T &f, typename FUNC_T::x_type a, typename FUNC_T::x_type b,
		                                         const typename FUNC_T::x_type &epsilon, PRINT_T p) {
    size_t k = 0; // search times
	typename FUNC_T::x_type x(a); 
	typename FUNC_T::y_type f_a = f(a), f_b = f(b);
	
	if (f_a * f_b >= 0)
		return Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type>(false, a);

	p << "k\ta\tb\tx\tf(a)\tf(b)\tf(x)\n";
	while (!((b - a) / 2 < epsilon)) {
		x = (a + b) / 2;
		typename FUNC_T::y_type f_x = f(x);
		p << k << "\t" << a << "\t" << b << "\t" << x << "\t" << f_a << "\t" << f_b << "\t" << f_x << '\n';
		if (f_a * f_x < 0) {
			b = x;
			f_b = f_x;
		}
		else { 
			a = x;
			f_a = f_x; 
		}
		++k;
	}
	return Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type>(true, (a + b) / 2);
}

} // namespace ends here

#endif

/* Test Case & User Manual

F should be f(x) whose root you are searching for.

class F {
public:
	typedef double x_type;
	typedef double y_type;

	double operator () (double x) const {
	    return x * x * x - x - 1;
	}
};

int main() {
	Felix::Non_linear_equations_search_returnT<double> res = 
		Felix::bisection(Felix::make_function(F()), 1.0, 2.0, 0.001, Felix::Stream_Printer(std::cout));
	if (res.get_valid()) {
		std::cout << "Answer: " << res.get_x() << std::endl;
	}
	else {
		std::cerr << "Root not found." << std::endl;
	}

	return 0;
}

*/

