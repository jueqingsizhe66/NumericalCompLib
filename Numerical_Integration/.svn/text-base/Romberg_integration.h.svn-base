#ifndef __ROMBERG_INTEGRATION_H__
#define __ROMBERG_INTEGRATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Romberg builds a series of sequences with increasing speed to converge,
// which can compute integration based on epsilon, the precision provided
// Algorithm:
// (1) N = 1, given epsilon, h_1 = b - a;
// (2) Compute:
//     T_1(0) = h_1 / 2 * [f(a) + f(b)]
// (3) h_2N = 1/2 * h_N
//     Compute:
//     T_2N(0) = 1/2 * T_N(0) + h_2N * SUM(1<=k<=N)(f(a + (2*k - 1)*h_2N))
// (4) M=N; N=2*N; k=1
// (5) Compute:
//     T_M(k) = (4^k * T_2M(k-1) - T_M(k-1)) / (4^k - 1)
// (6) if (M==1) goto (7) else M=M/2; ++k; goto (5)
// (7) if (|T_1(k) - T_1(k-1)| < epsilon) stop; return T_1(k) else goto (3)

/* Includes */
#include <vector>
using std::vector;

#include "../Basic_Components/basic_widgets.h"
#include "../Basic_Components/function.h"

/* Declaration */
namespace Felix { // namespace begins here

// h_2N = 1/2 * h_N
// T_2N = 1/2 * T_N + h_2N * SUM(1<=k<=N)(f(a+(2*k-1)*h_2N))
template <typename FUNC_T>
typename FUNC_T::y_type compute_T_2N(const FUNC_T &f, const typename FUNC_T::y_type &T_N, size_t N, typename FUNC_T::x_type h_N,
									 const typename FUNC_T::x_type &a) {
	h_N *= typename FUNC_T::x_type(0.5);
	typename FUNC_T::x_type sum(0);

	for (size_t k = 1; k <= N; ++k) {
	    sum += f(a + typename FUNC_T::x_type (2 * k - 1) * h_N);
	}
	return typename FUNC_T::y_type(0.5) * T_N + typename FUNC_T::x_type(h_N * sum);
}

// T_M(k) = (4^k * T_2M(k-1) - T_M(k-1)) / (4^k - 1)
template <typename T>
T compute_T_M_new(const T &T_2M, const T &T_M, int k) {
    T k_exp = pow(T(4), k);
	return (k_exp * T_2M - T_M) / (k_exp - 1);
}


template <typename FUNC_T>
typename FUNC_T::y_type Romberg_integration(const FUNC_T &f, const typename FUNC_T::x_type &a, const typename FUNC_T::x_type &b,
											const typename FUNC_T::x_type &epsilon) {
    size_t N = 1;
	int k = 0;
	typename FUNC_T::x_type h_N = b - a;
    
	vector<typename FUNC_T::y_type> T_N_seq, T_2N_seq;
	// compute T_1(0) and T_2(0)
	T_N_seq.push_back(h_N / typename FUNC_T::x_type(2) * (f(a) + f(b))); 
	T_2N_seq.push_back(compute_T_2N(f, T_N_seq[0], N, h_N, a));
	N *= 2;
	h_N *= typename FUNC_T::x_type(0.5);

	// compute T_1(1)
	k = 1;
	T_2N_seq.push_back(compute_T_M_new(T_2N_seq[0], T_N_seq[0], k));
	T_N_seq.push_back(typename FUNC_T::y_type(0));

	int limit = 1;
	while (Felix::abs(T_2N_seq[T_2N_seq.size() - 1] - T_N_seq[T_N_seq.size() - 2]) > epsilon) {
		// compute T_2N(0)
		typename FUNC_T::y_type tmp_T_2N = compute_T_2N(f, T_2N_seq[0], N, h_N, a);
		N *= 2;
		h_N *= typename FUNC_T::x_type(0.5);

		// refresh
		T_N_seq[0] = T_2N_seq[0];
		T_2N_seq[0] = tmp_T_2N;

		for (k = 1; k <= limit; ++k) {
		    typename FUNC_T::y_type tmp_T_N_new = compute_T_M_new(T_2N_seq[k - 1], T_N_seq[k - 1], k);
			T_N_seq[k] = T_2N_seq[k];
			T_2N_seq[k] = tmp_T_N_new; 
		}

		// deal with the last one
		T_2N_seq.push_back(compute_T_M_new(T_2N_seq[k - 1], T_N_seq[k - 1], k));
	    T_N_seq.push_back(typename FUNC_T::y_type(0));

		++limit;
	}
    return T_2N_seq[T_N_seq.size() - 1];
}


} // namespace ends here

#endif

/* Test Case & User Manual

class F {
public:
    typedef double x_type;
	typedef double y_type;
	double operator()(double x) const {
		return 4.0 / (1.0 + x * x);
	}
};

int main() {	
	std::cout << Felix::Romberg_integration(Felix::make_function(F()), 0.0, 1.0, 0.00001);

	return 0;
}

*/

