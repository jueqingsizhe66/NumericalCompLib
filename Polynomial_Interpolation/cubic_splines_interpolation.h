#ifndef __CUBIC_SPLINES_INTERPOLATION_H__
#define __CUBIC_SPLINES_INTERPOLATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Cubic splines interpolation is one of the most important interpolations
// in numerical analysis. 
// Given N+1 x and f(x) pair, S(x) suffices:
// (1) S(x) is a polynomial less than 3 times at every interval [x_k-1,x_k]
// (2) S(x_k) = y_k, 0<=k<=n
// (3) S(x) has continuous 2 derivative at [a=x_0,b=x_n]
// Algorithm:
// check any textbook about numerical analysis
// here, the implementation is about the first boundary condition: given S'(x_0) and S'(x_n)

/* Includes */
#include "../Basic_Components/matrix.h" 
#include "../Linear_Systems/Iterative_Methods/Gauss_Seidel_iteration.h"

#include <utility>
using std::pair;
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace begins here

// return type
template <typename T>
class cubic_splines_interpolation_polynomial {
public:
	cubic_splines_interpolation_polynomial(const vector<T> &x, const vector<T> &y, const vector<T> &h, const vector<T> &M) :
	    x_sequence_(x), y_sequence_(y), h_sequence_(h), M_sequence_(M) { }

    T operator () (const T &x) const {
        return this->get_y(x);																		
    }

	T get_y(const T &x) const {
		// find interval which x belongs to
		size_t lb = 0;
		bool find_lb = false;
		for (size_t i = 0; i < this->x_sequence_.size() - 1; ++i) {
			if (this->x_sequence_[i] <= x && x <= this->x_sequence_[i + 1]) {
			    lb = i;
				find_lb = true;
			}
		}
		if (!find_lb) {
			throw computation_exception("CUBIC SPLINES INTERPOLATION X IS NOT IN ANY INTERVAL!");
		}
		return this->M_sequence_[lb] / (T(6) * this->h_sequence_[lb]) * 
			   (this->x_sequence_[lb + 1] - x) * (this->x_sequence_[lb + 1] - x) * (this->x_sequence_[lb + 1] - x)
			   +
			   this->M_sequence_[lb + 1] / (T(6) * this->h_sequence_[lb]) *
			   (x - this->x_sequence_[lb]) * (x - this->x_sequence_[lb]) * (x - this->x_sequence_[lb])
			   +
			   T(1) / this->h_sequence_[lb] * 
			   (this->y_sequence_[lb + 1] - (this->M_sequence_[lb + 1] * this->h_sequence_[lb] * this->h_sequence_[lb] / T(6))) *
			   (x - this->x_sequence_[lb])
			   +
			   T(1) / this->h_sequence_[lb] *
			   (this->y_sequence_[lb] - (this->M_sequence_[lb] * this->h_sequence_[lb] * this->h_sequence_[lb] / T(6))) *
			   (this->x_sequence_[lb + 1] - x);
	}
public:
	vector<T> x_sequence_;
	vector<T> y_sequence_;
	vector<T> h_sequence_;
	vector<T> M_sequence_;
};


// Cubic splines interpolation
// x_y_vec: x and f(x) pairs
// dy_p: S'(x_0) and S'(x_n) pair
template <typename T>
cubic_splines_interpolation_polynomial<T> cubic_splines_interpolation(const vector<pair<T, T> > &x_y_vec, const pair<T, T> &dy_p) {
	matrix<T> tri_matrix(x_y_vec.size(), x_y_vec.size()),
		      d_matrix(x_y_vec.size(), 1);
	vector<T> x_seq, y_seq, h_seq;

	for (size_t i = 0; i != x_y_vec.size(); ++i) {
	    x_seq.push_back(x_y_vec[i].first);
		y_seq.push_back(x_y_vec[i].second);
		if (i > 0) {
		    h_seq.push_back(x_y_vec[i].first - x_y_vec[i - 1].first);
		}
	}

	// build tri_matrix and d_matrix
	for (size_t i = 0; i != x_y_vec.size(); ++i) {
		if (i < x_y_vec.size() - 1) {
			// set lower-left diagonal
			if (i == x_y_vec.size() - 2) {
			    tri_matrix[i + 1][i] = T(1);
			}
			else {
			    tri_matrix[i + 1][i] = h_seq[i] / (h_seq[i] + h_seq[i + 1]);
			}

			// set upper-right diagonal
			if (i == 0) {
			    tri_matrix[0][1] = T(1);
			}
			else {
			    tri_matrix[i][i + 1] = T(1) - tri_matrix[i][i - 1];
			}
		}
		tri_matrix[i][i] = T(2);
		if (i == 0) {
		    d_matrix[0][0] = T(6) / h_seq[0] * (((y_seq[1] - y_seq[0]) / h_seq[0]) - dy_p.first);
		}
		else if (i == x_y_vec.size() - 1) {
		    d_matrix[i][0] = T(6) / h_seq[i - 1] * (dy_p.second - ((y_seq[i] - y_seq[i - 1]) / h_seq[i - 1]));
		}
		else {
		    d_matrix[i][0] = T(6) / (h_seq[i - 1] + h_seq[i]) * ((y_seq[i + 1] - y_seq[i]) / h_seq[i] - (y_seq[i] - y_seq[i - 1]) / h_seq[i - 1]);
		}
	}

	// tri_matrix * M = d_matrix, solve M
	Jacobi_Gauss_Seidel_iteration_returnT<matrix<T> > res = Gauss_Seidel_iteration(tri_matrix, d_matrix, matrix<T>::epsilon, RAND_MAX, Null_Printer());
	vector<T> M;
	for (size_t i = 0; i != res.get_x().get_row_size(); ++i)
		M.push_back(res.get_x()[i][0]);
    return cubic_splines_interpolation_polynomial<T>(x_seq, y_seq, h_seq, M);
}


} // namespace ends here

#endif

/* Test Case & User Manual

int main() {
	double x_arr[4] = { 1, 4, 9, 16 };
	double y_arr[4] = { 1, 2, 3, 4 };
	std::pair<double, double> dy_p(0.5, 0.125);
	std::vector<std::pair<double, double> > x_y;

	for (int i = 0; i != 4; ++i)
		x_y.push_back(std::make_pair(x_arr[i], y_arr[i]));

	Felix::cubic_splines_interpolation_polynomial<double> res = Felix::cubic_splines_interpolation(x_y, dy_p);
	std::cout << res.get_y(5);

	return 0;
}

*/

