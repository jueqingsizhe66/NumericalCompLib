#ifndef __CONST_NUMBERS_AND_METRICS_H__
#define __CONST_NUMBERS_AND_METRICS_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// This class provides a set of const numbers that play 
// important roles in math system, such e and pi

/* Includes */
#include <string>
using std::string;

/* Declaration */
namespace Felix { // namesapce begins here

template <typename T>
class const_numbers_and_metrics {
public:
	static const T           e;
	static const T           pi;
	static T                 epsilon;
};

// static member data
template <typename T>
const T const_numbers_and_metrics<T>::e  = 2.7182818;

template <typename T>
const T const_numbers_and_metrics<T>::pi = 3.1415927;

template <typename T>
T const_numbers_and_metrics<T>::epsilon = 0.00001;

// template specializations
// int
template <>
const int const_numbers_and_metrics<int>::e  = 3;

template <>
const int const_numbers_and_metrics<int>::pi = 3;

template <>
int const_numbers_and_metrics<int>::epsilon = 1;


// unsigned int
template <>
const unsigned int const_numbers_and_metrics<unsigned int>::e  = 3;

template <>
const unsigned int const_numbers_and_metrics<unsigned int>::pi = 3;

template <>
unsigned int const_numbers_and_metrics<unsigned int>::epsilon = 1;


// std::string
template <>
const string const_numbers_and_metrics<string>::e  = "2.7182818";

template <>
const string const_numbers_and_metrics<string>::pi = "3.1415927";

template <>
string const_numbers_and_metrics<string>::epsilon = "0.00001";


} // namespace ends here

#endif

