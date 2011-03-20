#ifndef __ITERATIVE_SEARCH_H__
#define __ITERATIVE_SEARCH_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Iterative method is an effective way to solve non-linear equation
// Algorithm:
// Given f(x), x = f(x)
// (1) Given x_0, N, epsilon, and f(x).  k = 0;
// (2) x_k+1 = f(x_k);
// (3) if (|x_k+1 - x_k| < epsilon) stop;
// (4) if (k == N) stop; else ++k;
// Moreover, it can print all the computing procedures

/* Includes */
#include "../Auxiliary_Components/Printer.h"

/* Declaration */
namespace Felix { // namespace begins here

// return type
template <typename X_T>
class Non_linear_equations_search_returnT {
public:
	Non_linear_equations_search_returnT(bool v, const X_T &x) : valid_(v), x_(x) { }
	bool get_valid() const { return this->valid_; }
	X_T &get_x() { return this->x_; }
private:
	bool valid_;
	X_T  x_; // root
};


template <typename FUNC_T, typename CNT_T, typename PRINT_T>
Felix::Non_linear_equations_search_returnT<typename FUNC_T::x_type> iterative_search(const FUNC_T &f, typename FUNC_T::x_type x_k, const CNT_T &N, 
										 const typename FUNC_T::x_type &epsilon, PRINT_T p) {
	CNT_T k(0);
	typename FUNC_T::x_type x_k_1(x_k);
	typename FUNC_T::x_type diff(epsilon + 1);

	p << "k\tx\n";
	while (k < N && (!(Felix::abs(diff) < epsilon))) {
		p << k << '\t' << x_k << '\n';
		x_k_1 = f(x_k);
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
		return pow(Felix::const_numbers_and_metrics<double>::e, -1 * x);
	}
};

int main() {
	Felix::Non_linear_equations_search_returnT<double> res = 
		Felix::iterative_search(Felix::make_function(F()), 0.5, 20, 0.0001, Felix::Stream_Printer(std::cout));
	if (res.get_valid()) {
		std::cout << "Answer: " << res.get_x() << std::endl;
	}
	else {
		std::cerr << "Root not found." << std::endl;
	}

	return 0;
}

*/

