#ifndef __EULER_METHOD_H__
#define __EULER_METHOD_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Such equations as
// dy/dx = f(x, y)
// y(x_0) = y_0
// are called constant differential equations 
// Euler's method is to solve this kind of problem, 
// by computing the sequence of y(x) for different x
// Algorithm:
// here, an improved version of Euler's method is discussed
// (1) y_n+1(E) = y_n + h*f(x_n, y_n)
// (2) y_n+1 = y_n + h/2 * [f(x_n, y_n) + f(x_n+1, y_n+1(E))]

/* Includes */
#include <iostream>
using std::ostream;
#include <utility>
using std::pair;
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace begins here

// forward declaration
template <typename X_T, typename Y_T>
class Differential_equations_returnT;

template <typename X_T, typename Y_T>
ostream &operator << (ostream &, const Differential_equations_returnT<X_T, Y_T> &);


// return type
template <typename X_T, typename Y_T = X_T>
class Differential_equations_returnT {
public:
	Differential_equations_returnT(const vector<pair<X_T, Y_T> > &x_y_seq, const X_T &step) : x_y_seq_(x_y_seq), step_(step) { }

	size_t get_sequence_size() const { return this->x_y_seq_.size(); }
	Y_T get_y(const X_T &x) const {
		for (size_t i = 0; i != this->x_y_seq_.size(); ++i) {
		    if (this->x_y_seq_[i].first == x)
				return this->x_y_seq_[i].second;
		}
		throw computation_exception("EULER'S METHOD X-Y NOT FOUND!");
	}
	vector<pair<X_T, Y_T> > &get_x_y_sequence() { return this->x_y_seq_; }

	X_T &get_step() { return this->step_; }

	// ostream operation
	friend ostream &operator << <X_T, Y_T> (ostream &os, const Differential_equations_returnT<X_T, Y_T> &v);
private:
	vector<pair<X_T, Y_T> > x_y_seq_;
	X_T step_;
};

// implementation
template <typename X_T, typename Y_T>
ostream &operator << (ostream &os, const Differential_equations_returnT<X_T, Y_T> &v) {
	os << "step: " << v.step_ << '\n';
	for (size_t i = 0; i != v.x_y_seq_.size(); ++i) {
	    os << i << '\t' << v.x_y_seq_[i].first << '\t' << v.x_y_seq_[i].second << '\n'; 
	}
	return os;
}


// FUNC_T f is f(x, y)
template <typename FUNC_T>
Differential_equations_returnT<typename FUNC_T::x_type::first_type, typename FUNC_T::x_type::second_type> 
    Euler_method(const FUNC_T &f, const typename FUNC_T::x_type::first_type &x_0, 
													  const typename FUNC_T::x_type::second_type &y_0, 
													  const typename FUNC_T::x_type::first_type &x_n, 
													  const typename FUNC_T::x_type::first_type &step) {
    vector<pair<typename FUNC_T::x_type::first_type, typename FUNC_T::x_type::second_type> > seq;
	typename FUNC_T::x_type::second_type tmp_y, y_i;

	seq.push_back(std::make_pair(x_0, y_0));

	for (typename FUNC_T::x_type::first_type x_i = x_0 + step; x_i <= x_n; x_i += step) {
		tmp_y = seq[seq.size() - 1].second + step * f(std::make_pair(seq[seq.size() - 1].first, seq[seq.size() - 1].second));
		y_i = seq[seq.size() - 1].second + 
			  step * typename FUNC_T::x_type::first_type(0.5) * (f(std::make_pair(seq[seq.size() - 1].first, seq[seq.size() - 1].second)) + 
			                                                     f(std::make_pair(x_i, tmp_y)));
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
	Felix::Differential_equations_returnT<double> res = Felix::Euler_method(Felix::make_function(F()), 0.0, 1.0, 0.1, 0.02);
	std::cout << "step: " << res.get_step() << '\n'
	          << "sequence size:  " << res.get_sequence_size() << '\n';
    
	std::cout << "sequence: \n";
    for (unsigned int i = 0; i != res.get_sequence_size(); ++i)
		std::cout << res.get_x_y_sequence()[i].first << "\t" << res.get_x_y_sequence()[i].second << '\n';

	std::cout << "y(0.1) = " << res.get_y(0.1);

	return 0;
}

*/

