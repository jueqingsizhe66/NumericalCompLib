#ifndef __HERMITE_INTERPOLATION_H__
#define __HERMITE_INTERPOLATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Hermite interpolation can build a (2N+1)-polynomial H_2n+1 from N+1
// different pairs: y_k = f(x_k), 0<=k<=n
// and N+1 different pairs: y_k' = f'(x_k)
// where H_2n+1(x_k) = y_k and H_2n+1'(x_k) = y_k', 0<=k<=n 
// Algorithm:
// H_2n+1(x) = SUM(0<=k<=n){ y_k*(1-2*l_k'(x_k)*(x-x_k))*l_k(x)^2 + y_k'*(x-x_k)*l_k(x)^2 }
// l_k(x) = (x-x_0)*...*(x-x_k-1)*(x-x_k+1)*...*(x-x_n) / (x_k-x_0)*...*(x_k-x_k-1)*(x_k-x_k+1)*...*(x_k-x_n)
// l_k'(x) = SUM(0<=i<=n & i != k) { 1/(x_k-x_i) }

/* Includes */
#include <utility>
using std::pair;
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace begins here

// return type
template <typename T>
class Hermite_interpolation_polynomial {
public:
	Hermite_interpolation_polynomial(const vector<T> &x, const vector<T> &y, const vector<T> &dy, const vector<T> &l, 
		const vector<T> &dl) : x_sequence_(x), y_sequence_(y), dy_sequence_(dy), l_const_sqr_sequence_(l), dl_sequence_(dl) { }

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
			T tmp = prod * prod * this->l_const_sqr_sequence_[k];
		    sum += this->y_sequence_[k] * (T(1) - T(2) * this->dl_sequence_[k] * (x - this->x_sequence_[k])) * tmp 
				   + this->dy_sequence_[k] * (x - this->x_sequence_[k]) * tmp; 
		}
		return sum;
	}
public:
	vector<T> x_sequence_;
	vector<T> y_sequence_;
	vector<T> dy_sequence_;
	vector<T> l_const_sqr_sequence_;
	vector<T> dl_sequence_;
};


// Hermite interpolation
// x_y_dy_vec: pair<T, pair<T, T> >, they are x, y and y'
template <typename T>
Hermite_interpolation_polynomial<T> Hermite_interpolation(const vector<pair<T, pair<T, T> > > &x_y_dy_vec) {
	vector<T> x_seq, y_seq, dy_seq, l_const_sqr_seq, dl_seq;
	for (size_t k = 0; k != x_y_dy_vec.size(); ++k) {
	    x_seq.push_back(x_y_dy_vec[k].first);
		y_seq.push_back(x_y_dy_vec[k].second.first);
		dy_seq.push_back(x_y_dy_vec[k].second.second);

		T x_diff_prod(1), x_diff_sum(0);
		for (size_t j = 0; j != x_y_dy_vec.size(); ++j) {
			if (j != k) {
			    x_diff_prod *= x_y_dy_vec[k].first - x_y_dy_vec[j].first;
				x_diff_sum += T(1) / (x_y_dy_vec[k].first - x_y_dy_vec[j].first);
			}
		}
		l_const_sqr_seq.push_back((T(1) / x_diff_prod) * (T(1) / x_diff_prod));
		dl_seq.push_back(x_diff_sum);
	}
	return Hermite_interpolation_polynomial<T>(x_seq, y_seq, dy_seq, l_const_sqr_seq, dl_seq);
}


} // namespace ends here

#endif

/* Test Case & User Manual

int main() {
	double x_arr[2] = { 4, 9 };
	double y_arr[2] = { 2, 3 };
	double dy_arr[2] = { 0.25, 1.0/6.0 };
	std::vector<std::pair<double, std::pair<double, double> > > x_y_dy;

	for (int i = 0; i != 2; ++i)
		x_y_dy.push_back(std::make_pair(x_arr[i], std::make_pair(y_arr[i], dy_arr[i])));

	Felix::Hermite_interpolation_polynomial<double> res = Felix::Hermite_interpolation(x_y_dy);
	std::cout << res.get_y(5);

	return 0;
}

*/

