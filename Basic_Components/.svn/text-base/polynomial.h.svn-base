#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Polynomial is a basic form of function
// P(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n

/* Includes */
#include <vector>
using std::vector;

/* Declaration */
namespace Felix { // namespace Felix begins here

// forward declaration
template <typename T>
class polynomial;

template <typename T>
ostream &operator << (ostream &os, const polynomial<T> &p);

template <typename T>
istream &operator >> (istream &is, polynomial<T> &p);

template <typename T>
class polynomial {
public:
	polynomial() { }
	polynomial(const vector<T> &c) : coefficient_(c) { }

	T operator () (const T &x) const {
	    return this->compute(x);
	}

	T compute(const T &x) const {
		T sum(0), exp_x(1);
		for (size_t i = 0; i != this->coefficient_.size(); ++i) {
		    sum += this->coefficient_[i] * exp_x;
			exp_x *= x;
		}
		return sum;
	}

	// iostream operations
	friend ostream &operator << <T> (ostream &os, const polynomial<T> &p);
	friend istream &operator >> <T> (istream &is, polynomial<T> &p);
private:
    vector<T> coefficient_;
};


// IO Operations
template <typename T>
ostream &operator << (ostream &os, const polynomial<T> &p) {
	for (size_t i = 0; i != p.coefficient_.size(); ++i) {
		os << p.coefficient_[i] << "x^" << i << ((i == p.coefficient_.size() - 1) ? "" : " + ");
	}
	return os;
}

template <typename T>
istream &operator >> (istream &is, polynomial<T> &p) {
	T coeff;
	p.coefficient_.clear();
	while (is >> coeff) {
	    p.coefficient_.push_back(coeff);
	}
	return is;
}


} // namespace ends here

#endif

