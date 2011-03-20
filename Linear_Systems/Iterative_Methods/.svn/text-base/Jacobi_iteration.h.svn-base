#ifndef __JACOBI_ITERATION_H__
#define __JACOBI_ITERATION_H__

/* Basic Information */
// Author: Felix Huang
// School: BJUT

/* Description */
// Jacobi iteration is a simple method to solve large scale equation groups
// using iterative approach
// Algorithm:
// (1) Given A, b(Ax = b) and x_0; k = 0, epsilon and N
// (2) for 1<=i<=n  x_i(k+1) = (1 / a_ii) (b_i - SUM(1<=j<=n && j!=i)(a_ij * x_j(k)))
// (3) if (||x(k+1) - x(k)||_infinite < epsilon) stop; return x(k+1);
// (4) if (k==N) stop; else ++k;

/* Includes */


/* Declaration */
namespace Felix { // namespace begins here

// return type for both Jacobi iteration and Gauss-Seidel iteration
template <typename A_T, typename b_T = A_T>
class Jacobi_Gauss_Seidel_iteration_returnT {
public:
	Jacobi_Gauss_Seidel_iteration_returnT(bool v, const b_T &x) :
	  valid_(v), x_(x) { }
	bool                   &get_valid() { return this->valid_; }
	b_T                    &get_x() { return this->x_; }
private:
    bool                   valid_;
	b_T                    x_; // the root of the equations
};


// implementation of step iteration, suitable for both Jacobi iteration and Gauss-Seidel iteration
template <typename A_T>
typename A_T::ele_type Jacobi_Gauss_Seidel_step_iteration(const A_T &a, const A_T &b, const A_T &x_k, A_T &x_k_1, const typename A_T::ele_type &w) {
	if ((x_k.get_row_size() != x_k_1.get_row_size()) || (x_k.get_col_size() != x_k_1.get_col_size()))
		throw "JACOBI GAUSS SEIDEL STEP ITERATION X SIZE UNMATCHED!";

	A_T old_x_k(x_k);
	for (size_t i = 0; i != a.get_row_size(); ++i) {
		typename A_T::ele_type sum(0);
		for (size_t j = 0; j != a.get_row_size(); ++j) {
		    if (j != i)
				sum += a[i][j] * x_k[j][0];
		}
		x_k_1[i][0] = (typename A_T::ele_type(1) - w) * x_k[i][0] + w * (b[i][0] - sum) / a[i][i];
	}
	return (x_k_1 - old_x_k).get_infinite_norm();
}


// Jacobi iteration 
template <typename A_T, typename CNT_T, typename PRINT_T>
Jacobi_Gauss_Seidel_iteration_returnT<A_T> Jacobi_iteration(const A_T &a, const A_T &b, const typename A_T::ele_type &epsilon, const CNT_T &N, PRINT_T p) {
    CNT_T k(0);
	A_T x_k(a.get_row_size(), 1);
	A_T x_k_1(x_k); 
	typename A_T::ele_type diff(epsilon + typename A_T::ele_type(1));

	while (k < N && diff > epsilon) {
		p << "k:  " << k << '\n';
		p << "x:\n" << x_k;
		(k == 0) ? (p << "diff:  ---\n\n") : (p << "diff:  " << diff << "\n\n");

		// use different x sequence
		diff = Felix::Jacobi_Gauss_Seidel_step_iteration(a, b, x_k, x_k_1, typename A_T::ele_type(1));
		x_k = x_k_1;
		++k;
	} 
	if (k == N)
		return Felix::Jacobi_Gauss_Seidel_iteration_returnT<A_T>(false, x_k_1);
	else
		return Felix::Jacobi_Gauss_Seidel_iteration_returnT<A_T>(true, x_k_1);
}


} // namespace ends here


#endif

/* Test Case & User Manual

int main() {
	double A_arr[3][3] = { {11, -3,  -2},
	                       {-1,  5,  -3},
	                       {-2, -12, 19} };
	double b_arr[3][1] = { {3},
		                   {6},
						   {-7} };
	Felix::matrix<double> A(A_arr), b(b_arr); 

	Felix::Jacobi_Gauss_Seidel_iteration_returnT<Felix::matrix<double> > res = 
		Felix::Jacobi_iteration(A, b, 0.001, 30, Felix::Stream_Printer(std::cout));

	if (res.get_valid()) {
		std::cout << res.get_x();
	}
	else {
		std::cerr << "Root not found.";
	}

	return 0;
}

*/

