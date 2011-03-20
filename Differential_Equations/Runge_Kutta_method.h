#ifndef __RUNGE_KUTTA_METHOD_H__
#define __RUNGE_KUTTA_METHOD_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Runge Kutta method is an efficient way to solve constant 
// differential equations
// Algorithm:
// y_n+1 = y_n + 1/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
// k_1 = h * f(x_n, y_n)
// k_2 = h * f(x_n + 1/2 * h, y_n + 1/2 * k_1)
// k_3 = h * f(x_n + 1/2 * h, y_n + 1/2 * k_2)
// k_4 = h * f(x_n + h, y_n + k_3)

/* Includes */
#include "Euler_method.h"

/* Declaration */
namespace Felix { // namespace begins here

// FUNC_T f is f(x, y)
template <typename FUNC_T>
Differential_equations_returnT<typename FUNC_T::x_type::first_type, typename FUNC_T::x_type::second_type> 
    Runge_Kutta_method(const FUNC_T &f, const typename FUNC_T::x_type::first_type &x_0, 
													  const typename FUNC_T::x_type::second_type &y_0, 
													  const typename FUNC_T::x_type::first_type &x_n, 
													  const typename FUNC_T::x_type::first_type &step) {
    vector<pair<typename FUNC_T::x_type::first_type, typename FUNC_T::x_type::second_type> > seq;
	typename FUNC_T::x_type::second_type y_i, k1, k2, k3, k4;

	seq.push_back(std::make_pair(x_0, y_0));

	for (typename FUNC_T::x_type::first_type x_i = x_0 + step; x_i <= x_n; x_i += step) {
		k1 = step * f(std::make_pair(seq[seq.size() - 1].first, seq[seq.size() - 1].second));
		k2 = step * f(std::make_pair(seq[seq.size() - 1].first + step * typename FUNC_T::x_type::first_type(0.5), 
			          seq[seq.size() - 1].second + k1 * typename FUNC_T::x_type::first_type(0.5)));
		k3 = step * f(std::make_pair(seq[seq.size() - 1].first + step * typename FUNC_T::x_type::first_type(0.5), 
			          seq[seq.size() - 1].second + k2 * typename FUNC_T::x_type::first_type(0.5)));
		k4 = step * f(std::make_pair(seq[seq.size() - 1].first + step, seq[seq.size() - 1].second + k3));
        y_i = seq[seq.size() - 1].second + typename FUNC_T::x_type::first_type(1) / typename FUNC_T::x_type::first_type(6) *
			(k1 + typename FUNC_T::x_type::first_type(2) * k2 + typename FUNC_T::x_type::first_type(2) * k3 + k4);

		seq.push_back(std::make_pair(x_i, y_i));
	}
	return Differential_equations_returnT<typename FUNC_T::x_type::first_type, typename FUNC_T::x_type::second_type>(seq, step);
}


} // namespace ends here

#endif

/* Test Case & User Manual

class F {
public:
	typedef std::pair<double, double> x_type;
	typedef double y_type;
	double operator()(x_type x) const {
		return x.first + x.second;
	}
};

int main() {
	Felix::Differential_equations_returnT<double> res = Felix::Runge_Kutta_method(Felix::make_function(F()), 0.0, 1.0, 0.1, 0.1);
	std::cout << "step: " << res.get_step() << '\n'
	          << "sequence size:  " << res.get_sequence_size() << '\n';
    
	std::cout << "sequence: \n";
    for (unsigned int i = 0; i != res.get_sequence_size(); ++i)
		std::cout << res.get_x_y_sequence()[i].first << "\t" << res.get_x_y_sequence()[i].second << '\n';

	std::cout << "y(0.1) = " << res.get_y(0.1);

	return 0;
}

*/

