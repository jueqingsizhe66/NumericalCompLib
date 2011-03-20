#ifndef __FUNCTION_H__
#define __FUNCTION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Function is an significant component of numerical algorithms
// Anything like y = f(x) is definded here, which
// can be used to build higher level system

namespace Felix { // namespace Felix begins here

template <typename FUNC_T>
class function {
public:
	// define types
	typedef FUNC_T f_type;
	typedef typename FUNC_T::x_type x_type;
	typedef typename FUNC_T::y_type y_type;

	// constructor
	function(const FUNC_T &f) : f_(f) { }

	// copy constructor --- default
	// operator = --- default
	// destructor --- default

	// y = f(x)
	y_type operator () (const x_type &x) const { return this->f_(x); } 
private:
    FUNC_T f_;
};

// convenient function
template <typename FUNC_T>
function<FUNC_T> make_function(const FUNC_T &f) { return function<FUNC_T>(f); }

} // namespace ends here

#endif

