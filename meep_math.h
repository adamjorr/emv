#ifndef __MEEP_MATH_INCLUDED__
#define __MEEP_MATH_INCLUDED__
#include <array>

namespace meep_math{
	double binomial_cdf(int successes, int trials, long double p);
	double binomial_pdf(int successes, int trials, long double p);
	int binomial_coeff(int n, int k);
	std::array<double,2> solve_quadratic(double a, double b, double c);
	std::array<double,4> solve_quartic(double a, double b, double c, double d, double e);
}
#endif