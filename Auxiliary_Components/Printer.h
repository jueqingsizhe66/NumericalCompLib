#ifndef __PRINTER_H__
#define __PRINTER_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Printers are used as components in algorithms to print
// the procedural results.  A variety of printers are provided
// to cater to different requirements

/* Includes */
#include <iostream>
using std::ostream;

/* Declaration */
namespace Felix { // namespace Felix begins here 

// print to stream
class Stream_Printer {
public:
	Stream_Printer(ostream &os = std::cout) : os_(os) { }
	Stream_Printer(const Stream_Printer &rhs) : os_(rhs.os_) { }

	template <typename T>
	Stream_Printer &operator << (const T &val) {
		this->os_ << val;
		return *this;
	}
private:
	ostream &os_;
};

// do nothing
class Null_Printer {
public:
	template <typename T>
	Null_Printer &operator << (const T &val) {
		return *this;
	}
};

// print std containers, facilitating debugging
class StdCon_Printer {
public:
	StdCon_Printer(ostream &os = std::cout) : p_(os) { }

    template <typename CON_T>
	StdCon_Printer &operator << (const CON_T &con) {
		this->p_ << "[ ";
		for (size_t i = 0; i != con.size(); ++i) {
			this->p_ << con[i] << (i == con.size() - 1 ? " " : ", ");
		}
		this->p_ << "]\n";
		return *this;
	}
private:
	Stream_Printer p_;
};

} // namespace ends here

#endif

