#ifndef __FUNCTION_ANALYZER_H__
#define __FUNCTION_ANALYZER_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// The analyzer is to parse an string input and build
// an inner type that represents the function embodied
// by the string

/* Includes */
#include <string>
using std::string;

#include "../Basic_Components/function.h"

/* Declaration */
namespace Felix { // namespace Felix begins here 

// NOTICE: This function has not been implemented yet.
// A dummy function, y = x, is used as return.
template <typename T>
class Dummy_function {
public:
	typedef T x_type;
	typedef T y_type;
	T operator() (const T &x) const {
	    return x;
	}
};


class function_analyzer {
public:
	template <typename T>
	static Felix::function<T> get_function_from_string(const string &f, const T &type_hint) {
		return Felix::make_function(Dummy_function<T>());
	}
};


} // namespace ends here


#endif

