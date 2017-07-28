#ifndef __MEEP_MATH_INCLUDED__
#define __MEEP_MATH_INCLUDED__
#include <array>
#include <functional>
#include <complex>

namespace meep_math{
	double binomial_cdf(int successes, int trials, long double p);
	double binomial_pdf(int successes, int trials, long double p);
	int binomial_coeff(int n, int k);
	std::array<double,2> solve_quadratic(double a, double b, double c);
	std::array<std::complex<double>,4> solve_quartic(double a, double b, double c, double d, double e);
	double nr_root(std::function<double(double)> f, std::function<double(double)> f_prime, double init, double tolerance=.001, int maxiter=1000);
}
#endif