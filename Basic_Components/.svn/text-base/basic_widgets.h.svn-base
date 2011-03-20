#ifndef __BASIC_WIDGETS_H__
#define __BASIC_WIDGETS_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Here provides a variety of functions that are used frequently

/* Includes */
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <string>

/* Declaration */
namespace Felix { // namespace Felix begins here.

// compute |val|
template <typename T>
T abs(const T &val) {
	return val >= 0 ? val : -1 * val;
}

// notice the interval is [lb, ub)
template <typename T>
T positive_rand(const T &lb, const T &ub) {
	static bool set_seed = false;
	if (!set_seed) {
	    srand(static_cast<unsigned int>(time(NULL)));
		set_seed = true;
	}
	double ratio = static_cast<double>(rand()) / RAND_MAX;
	return static_cast<T>(lb + (ub - lb) * ratio);
}

// uniform exception
class computation_exception : public std::exception {
public:
	computation_exception(const std::string &msg) : msg_(msg) { }
	~computation_exception() throw() { }

	virtual const char* what() const throw() {
	    return msg_.c_str();
	}
private:
	std::string msg_;
};


} // namespace ends here


#endif

