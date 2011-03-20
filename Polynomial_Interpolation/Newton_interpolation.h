#ifndef __NEWTON_INTERPOLATION_H__
#define __NEWTON_INTERPOLATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Newton interpolation can build a N-polynomial N_n from N+1
// different pairs: y_k = f(x_k), 0<=k<=n
// where L_n(x_k) = y_k, 0<=k<=n
// Algorithm:
// N_n(x) = f(x_0) + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + ...
// Mean Inequality:
// f[x_0,x] = (f(x)-f(x_0)) / (x-x_0)
// ...
// f[x_0,x_1,x_2,...,x_n] = (f[x_0,...,x_n-2,x_n]-f[x_0,...,x_n-2,x_n-1]) / (x_n-x_n-1)


/* Includes */
#include <utility>
using std::pair;
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace begins here

// return type
template <typename T>
class Newton_interpolation_polynomial {
public:
	Newton_interpolation_polynomial(const vector<T> &x, const vector<T> &ie_top, const vector<T> &ie_bot, const T &ie_n) :
	                               x_sequence_(x), mean_inequality_top_sequence_(ie_top), mean_inequality_bot_sequence_(ie_bot),
									   mean_inequality_n_(ie_n) { }

    T operator () (const T &x) const {
        return this->get_y(x);																		
    }

	T get_y(const T &x) const {
		T sum(0), x_prod(1);
		
		size_t i;
		for (i = 0; i != this->mean_inequality_top_sequence_.size(); ++i) {
		    sum += this->mean_inequality_top_sequence_[i] * x_prod;
			x_prod *= x - this->x_sequence_[i];
		}
		sum += this->mean_inequality_n_ * x_prod;
		return sum;
	}

	void add_x_y_pair(const pair<T, T> &x_y) {
		if (x_y.first <= this->x_sequence_[this->x_sequence_.size() - 1])
			throw computation_exception("NEWTON INTERPOLATION ADDING INVALID PAIR!");
		this->x_sequence_.push_back(x_y.first);
		this->mean_inequality_top_sequence_.push_back(this->mean_inequality_n_);
		this->mean_inequality_bot_sequence_.push_back(this->mean_inequality_n_);
		T ie_n(x_y.second);
		for (size_t i = 0; i != this->mean_inequality_bot_sequence_.size(); ++i) {
		    T tmp = (ie_n - this->mean_inequality_bot_sequence_[i]) / 
				    (this->x_sequence_[this->x_sequence_.size() - 1] - this->x_sequence_[this->x_sequence_.size() - 2 - i]);
			this->mean_inequality_bot_sequence_[i] = ie_n;
			ie_n = tmp;
		}
		this->mean_inequality_n_ = ie_n;
	}
public:
	vector<T> x_sequence_;
	vector<T> mean_inequality_top_sequence_;
	vector<T> mean_inequality_bot_sequence_;
	T mean_inequality_n_;
};


// Newton interpolation
// pair<T, T> first is x, second is y, namely f(x)
template <typename T>
Newton_interpolation_polynomial<T> Newton_interpolation(const vector<pair<T, T> > &x_y_vec) {
	vector<T> x_seq, ie_top, ie_bot, tmp;
	T ie_n(0);

	for (size_t i = 0; i != x_y_vec.size(); ++i) {
	    x_seq.push_back(x_y_vec[i].first);
		tmp.push_back(x_y_vec[i].second);
	}

	for (size_t i = x_seq.size() - 1; i >= 1; --i) {
		ie_top.push_back(tmp[0]);
		ie_bot.push_back(tmp[tmp.size() - 1]);
		for (size_t j = 1; j <= i; ++j) {
		    tmp[j - 1] = (tmp[j] - tmp[j - 1]) / (x_seq[j + x_seq.size() - 1 - i] - x_seq[j - 1]);
		}
		tmp.resize(i);
	}
	ie_n = tmp[0];

	return Newton_interpolation_polynomial<T>(x_seq, ie_top, ie_bot, ie_n);
}


} // namespace ends here

#endif

/* Test Case & User Manual

int main() {
	double x_arr[2] = { 1, 4 };
	double y_arr[2] = { 1, 2 };
	std::vector<std::pair<double, double> > x_y;

	for (int i = 0; i != 2; ++i)
		x_y.push_back(std::make_pair(x_arr[i], y_arr[i]));

	Felix::Newton_interpolation_polynomial<double> res = Felix::Newton_interpolation(x_y);
	std::cout << res.get_y(5) << std::endl;

	res.add_x_y_pair(std::make_pair(9, 3));
	std::cout << res.get_y(5);

	return 0;
}

*/

