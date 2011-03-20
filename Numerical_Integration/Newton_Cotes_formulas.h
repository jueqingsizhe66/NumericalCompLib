#ifndef __NEWTON_COTES_FORMULAS_H__
#define __NEWTON_COTES_FORMULAS_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
//  _b
// |   f(x) dx =(approximate) SUM(0<=k<=n)(A_k*f(x_k))
//_| a
//
// The essence of Newton-Cotes formulas is to substitute f(x)
// for L_n(x) where L_n(x) is the Lagrange interpolation polynomial
// for f(x). Specifically, in practice, n=4.
// Algorithm:
//
//  _b
// |   f(x) dx =(approximate) (b-a)/90 [ 7*f(a) + 32*f(a+h) + 12*f(a+2h) + 32*f(a+3h) + 7*f(b)]
//_| a
//
// where x_0=a, x_k=x_0 + k*h, 1<=k<=4; h=(b-a)/4

/* Includes */
#include "../Basic_Components/function.h"

/* Declaration */
namespace Felix { // namespace begins here

template <typename FUNC_T>
typename FUNC_T::y_type Newton_Cotes_quadrature_formula(const FUNC_T &f, const typename FUNC_T::x_type &a, const typename FUNC_T::x_type &b) {
	typename FUNC_T::x_type h = (b - a) / typename FUNC_T::x_type(4);
	return (b - a) / typename FUNC_T::x_type(90) * (
		    typename FUNC_T::x_type(7) * f(a) +
			typename FUNC_T::x_type(32) * f(a + h) +
			typename FUNC_T::x_type(12) * f(a + 2 * h) +
			typename FUNC_T::x_type(32) * f(a + 3 * h) +
			typename FUNC_T::x_type(7) * f(b)
		);
}

// Composite rules for Newton-Cotes quadrature formula
// h_n = (b-a)/N and in every interval, apply Newton-Cotes rule
template <typename FUNC_T, typename N_T>
typename FUNC_T::y_type Newton_Cotes_quadrature_formula(const FUNC_T &f, const typename FUNC_T::x_type &a, const typename FUNC_T::x_type &b,
														const N_T &N) {
	typename FUNC_T::x_type h_n = (b - a) / static_cast<typename FUNC_T::x_type> (N);
	typename FUNC_T::y_type sum1(0), sum2(0), sum3(0), sum4(0);
	
	for (N_T k = 0; k < N; ++k) {
	    sum1 += f(a + (N_T(4) * k + N_T(1)) * h_n / typename FUNC_T::x_type(4));
		sum2 += f(a + (N_T(4) * k + N_T(2)) * h_n / typename FUNC_T::x_type(4));
		sum3 += f(a + (N_T(4) * k + N_T(3)) * h_n / typename FUNC_T::x_type(4));
		sum4 += (k == N_T(0) ? typename FUNC_T::y_type(0) : f(a + (N_T(4) * k) * h_n / typename FUNC_T::x_type(4)));
	}

	return h_n / typename FUNC_T::x_type(90) * (
		    typename FUNC_T::x_type(7) * f(a) +
			typename FUNC_T::x_type(7) * f(b) +
			typename FUNC_T::x_type(32) * sum1 +
			typename FUNC_T::x_type(12) * sum2 +
			typename FUNC_T::x_type(32) * sum3 +
			typename FUNC_T::x_type(14) * sum4
		);
}


} // namespace ends here

#endif

/* Test Case & User Manual

class F {
public:
    typedef double x_type;
	typedef double y_type;
	double operator()(double x) const {
		return 1.0 / (1.0 + x * x);
	}
};

int main() {
	std::cout << Felix::Newton_Cotes_quadrature_formula(Felix::make_function(F()), 0.0, 1.0);

	return 0;
}

*/
