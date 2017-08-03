#include <string>
#include "meep_math.h"

namespace meep_math{
	//TODO: cache this somewhere or we're gonna be slow
	double binomial_cdf(int successes, int trials, long double p){
		long double cdf = 0.0;
		for(int i = 0; i < trials; i++){
			cdf += meep_math::binomial_pdf(i, trials, p);
		}
		return cdf;
	}

	double binomial_pdf(int successes, int trials, long double p){
		return (double)meep_math::binomial_coeff(trials, successes) * pow(p, successes) * pow(1 - p,trials - successes);
	}

	//n choose k
	int binomial_coeff(int n, int k){
		if (k > n){
			throw std::runtime_error("error calculating binomial with n = " + std::to_string(n) + " and k = " + std::to_string(k) );
		}
		int product = 1;
		for(int i=1; i <= k; i++){
			product *= ((n + 1 - i) / i);
		}
		return product;
	}

	std::array<double,2> solve_quadratic(double a, double b, double c){
		std::array<double,2> x = {0.0,0.0};
		x[0] = (-b - sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
		x[1] = (-b + sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
		return x;
	}

	//https://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	std::array<std::complex<double>,4> solve_quartic(double a, double b, double c, double d, double e){
		std::array<std::complex<double>,4> x = {0.0,0.0,0.0,0.0};
		double p = (8.0 * a * c - 3 * pow(b,2))/(8.0 * pow(a,2));
		double q = (pow(b,3) - 4.0 * a * b * c + 8.0 * pow(a,2) * d)/(8.0 * pow(a,3));
		double disc = 256 * pow(a,3) * pow(e,3) - 192 * pow(a,2) * b * d * pow(e,2) - 128 * pow(a,2) * pow(c,2) * pow(e,2)
						+ 144 * pow(a,2) * c * pow(d,2) * e - 27 * pow(a,2) * pow(d,4)
						+ 144 * a * pow(b,2) * c * pow(e,2) - 6 * a * pow(b,2) * pow(d,2) * e - 80 * a * b * pow(c,2) * d * e + 18 * a * b * c * pow(d,3)
						+ 16 * a * pow(c,4) * e - 4 * a * pow(c,3) * pow(d,2) - 27 * pow(b,4) * pow(e,2) + 18 * pow(b,3) * c * d * e - 4 * pow(b,3) * pow(d,3)
						- 4 * pow(b,2) * pow(c,3) * e + pow(b,2) * pow(c,2) * pow(d,2);
		double disc0 = pow(c,2) - 3.0 * b * d + 12.0 * a * e;
		double disc1 = 2.0 * pow(c,3) - 9 * b * c * d + 27 * pow(b,2) * e + 27 * a * pow(d,2) - 72 * a * c * e;
		std::complex<double> Q;
		if(disc != 0 && disc0 == 0){
			Q = std::pow(disc1,1/3); //special case
		}
		else{
			Q = std::pow((disc1 + std::sqrt(-27 * disc))/2,1/3);	
		}
		std::complex<double> S = 0.5 * std::sqrt(-(2.0/3.0)*p + 1.0/(3.0*a) * (Q + disc0/Q));

		//special case. when S == 0 we try other cube roots.
		//see https://en.wikipedia.org/wiki/Cube_root
		//if none of them work, it's because the depressed quartic is biquadratic
		if(S == 0.0){
			std::complex<double> old_Q = Q;
			Q = old_Q * std::complex<double>(-0.5,std::sqrt(3)/2);
			S = 0.5 * std::sqrt(-(2.0/3.0)*p + 1.0/(3.0*a) * (Q + disc0/Q));
			if (S == 0.0){
				Q = old_Q * std::complex<double>(-0.5,-std::sqrt(3)/2);
				S = 0.5 * std::sqrt(-(2.0/3.0)*p + 1.0/(3.0*a) * (Q + disc0/Q));
				if (S == 0.0){
					std::array<double,2> z = meep_math::solve_quadratic(a,c,e);
					x[0] = std::sqrt(z[0]);
					x[1] = -std::sqrt(z[0]);
					x[2] = std::sqrt(z[1]);
					x[3] = -std::sqrt(z[2]);
					return x;
				}
			}
		}
		//compute roots
		x[0] = -b / (4.0 * a) - S - 0.5 * std::sqrt(-4.0 * pow(S,2) - 2*p + q/S);
		x[1] = -b / (4.0 * a) - S + 0.5 * std::sqrt(-4.0 * pow(S,2) - 2*p + q/S);
		x[2] = -b / (4.0 * a) + S - 0.5 * std::sqrt(-4.0 * pow(S,2) - 2*p + q/S);
		x[3] = -b / (4.0 * a) + S + 0.5 * std::sqrt(-4.0 * pow(S,2) - 2*p + q/S);
		return x;
	}

	//newton-raphson root finding DO NOT USE, DOESN'T WORK
	double nr_root(std::function<double(double)> f, std::function<double(double)> f_prime, double init, double tolerance, int maxiter){
		double x = init;
		int iter = 0;
		double f_x = f(x);
		while (f_x > tolerance || iter++ <= maxiter){
			double f_pm = f_prime(x);
			if (f_pm == 0){
				throw std::runtime_error("Divison by zero; x = " + std::to_string(x) + " f(x) = " + std::to_string(f_x) + " f'(x) = " + std::to_string(f_pm));
			}
			x -= f_x / f_prime(x);
			f_x = f(x);
		}
		if (iter == maxiter){
			throw std::runtime_error("Newton-rafsom failed to converge");
		}
		return x;
	}



}