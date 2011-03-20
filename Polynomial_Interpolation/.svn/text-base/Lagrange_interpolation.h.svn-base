#ifndef __LAGRANGE_INTERPOLATION_H__
#define __LAGRANGE_INTERPOLATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Lagrange interpolation can build a N-polynomial L_n from N+1
// different pairs: y_k = f(x_k), 0<=k<=n
// where L_n(x_k) = y_k, 0<=k<=n
// Algorithm:
// L_n(x) = SUM(0<=k<=n)(y_k * PROD(0<=j<=n&&j!=k)((x - x_j)/(x_k - x_j)))
// L_n can be divided into two parts: one that is related to x, and the other
// that can be completely decided by input pairs
// In this implementation, the first, x-related part is called: x_sequence_
// the other is called y_and_const_sequence_

/* Includes */
#include <utility>
using std::pair;
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace begins here

// return type
template <typename T>
class Lagrange_interpolation_polynomial {
public:
	Lagrange_interpolation_polynomial(const vector<T> &x, const vector<T> &y_and_const) :
	                                                                    x_sequence_(x), y_and_const_sequence_(y_and_const) { }
    T operator () (const T &x) const {
        return this->get_y(x);																		
    }

	T get_y(const T &x) const {
		T sum(0);
		vector<T> x_minus_x_k_vec;

		for (size_t i = 0; i != this->x_sequence_.size(); ++i)
			x_minus_x_k_vec.push_back(x - this->x_sequence_[i]);

		for (size_t k = 0; k != this->x_sequence_.size(); ++k) {
			T prod(1);
			for (size_t j = 0; j != this->x_sequence_.size(); ++j) {
				if (j != k) {
				    prod *= x_minus_x_k_vec[j];
				}
			}
		    sum += (prod * this->y_and_const_sequence_[k]); 
		}
		return sum;
	}
public:
	vector<T> x_sequence_;
	vector<T> y_and_const_sequence_;
};


// Lagrange interpolation
// pair<T, T> first is x, second is y, namely f(x)
template <typename T>
Lagrange_interpolation_polynomial<T> Lagrange_interpolation(const vector<pair<T, T> > &x_y_vec) {
	vector<T> x_seq, y_and_const_seq;
	for (size_t k = 0; k != x_y_vec.size(); ++k) {
	    x_seq.push_back(x_y_vec[k].first);

		T x_diff_prod(1);
		for (size_t j = 0; j != x_y_vec.size(); ++j) {
			if (j != k) {
			    x_diff_prod *= x_y_vec[k].first - x_y_vec[j].first;
			}
		}
		y_and_const_seq.push_back(x_y_vec[k].second / x_diff_prod);
	}
	return Lagrange_interpolation_polynomial<T>(x_seq, y_and_const_seq);
}


} // namespace ends here

#endif

/* Test Case & User Manual

int main() {
	double x_arr[3] = { 1, 4, 9 };
	double y_arr[3] = { 1, 2, 3 };
	std::vector<std::pair<double, double> > x_y;

	for (int i = 0; i != 3; ++i)
		x_y.push_back(std::make_pair(x_arr[i], y_arr[i]));

	Felix::Lagrange_interpolation_polynomial<double> res = Felix::Lagrange_interpolation(x_y);
	std::cout << res.get_y(5);

	return 0;
}

*/

