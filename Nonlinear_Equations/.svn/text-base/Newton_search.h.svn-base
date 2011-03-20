#ifndef __NEWTON_SEARCH_H__
#define __NEWTON_SEARCH_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Newton method is named after one of the most preeminent mathematicians Isaac Newton
// It can sovle non-linear equation very effectively
// Algorithm:
// Given f(x), f(x) = 0
// (1) Given x_0, N, epsilon, f(x), f'(x).  k = 0;
// (2) x_k+1 = x_k - f(x_k) / f'(x_k);
// (3) if (|x_k+1 - x_k| < epsilon) stop;
// (4) if (k == N) stop; else { ++k; goto(2); }

/* Includes */
#include "iterative_search.h"
#include "../Basic_Components/Basic_Components.h"

/* Declaration */
namespace Felix { // namespace begins here

// Newton Function
template <typename FUNC_T, typename FUNC_D_T>
class Newton_Func {
public:
	typedef typename FUNC_T::x_type x_type;
	typedef typename FUNC_T::y_type y_type;
	Newton_Func(const FUNC_T &f, const FUNC_D_T &f_d) : f_(f), f_d_(f_d) { }
	x_type operator () (const x_type &x) const { return x - f_(x) / f_d_(x); }
private:
	FUNC_T f_;
	FUNC_D_T f_d_;
};

template <typename FUNC_T, typename FUNC_D_T, typename CNT_T, typename PRINT_T>
Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> Newton_search(const FUNC_T &f, const FUNC_D_T &f_derivatives, typename FUNC_T::x_type x_k, const CNT_T &N, 
	                                  const typename FUNC_T::x_type &epsilon, PRINT_T p) {
	return Felix::iterative_search(Felix::make_function(Newton_Func<FUNC_T, FUNC_D_T>(f, f_derivatives)), x_k, N, epsilon, p);
}


} // namesapce ends here


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
		Felix::Newton_search(Felix::make_function(F()), Felix::make_function(Fd()), 1.5, 20, 0.001, Felix::Stream_Printer(std::cout));
	if (res.get_valid()) {
		std::cout << "Answer: " << res.get_x() << std::endl;
	}
	else {
		std::cerr << "Root not found." << std::endl;
	}

	return 0;
}

*/

