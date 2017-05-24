#include "seqem.h"
#include "plpdata.h"
#include <string>
#include <iostream>
#include <cmath>

int main(int argc, char *argv[]){
	// Seqem seq("testdata/test.sam","testdata/test.fa");
	std::clog.precision(20);

	std::vector<char> x(9,'C');
	x.push_back('A');
	Pileupdata data(x);
	Seqem seq(data);

	// Genotype g("CC");
	// std::tuple<double> theta = std::make_tuple(0.01);
	// int counter = 0;
	// std::cout << "X = " << x << std::endl;
	// std::cout << "Theta = " << theta << "\t3 * Theta = " << std::get<0>(theta) * 3;
	// while (counter < 10){
	// 	double likelihood = 0.0;
	// 	std::vector<double> s(3,0.0);
	// 	for (Genotype g : Genotype::enumerate_gts(2)){
	// 		//calc likelihood
	// 		double b = seq.px_given_gtheta(x, g, theta); // log
	// 		double c = seq.pg(g); //log
	// 		likelihood += exp(b) * c;

	// 		//calc s
	// 		std::vector<double> site_s = seq.calc_s(x,g,theta);
	// 		double pg_x = seq.pg_x_given_theta(g,x,theta);
	// 		for(int i = 0; i < s.size(); ++i){
	// 			s[i] += pg_x * site_s[i];	
	// 		}
	// 	}

	// 	double a = (3.0 * s[0] + 3.0 * s[1] + s[2]);
	// 	double b = - (3.0/2 * s[0] + s[1] + 3.0/2 * s[2]);
	// 	double c = s[2] / 2;
	// 	double mu = (-b - sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
	// 	theta = std::make_tuple(mu);

	// 	std::cout << "\tLikelihood = " << likelihood << std::endl;
	// 	std::cout << "Theta = " << theta << "\t3 * Theta = " << std::get<0>(theta) * 3;
	// 	counter++;
	// }
	// std::cout << std::endl;
	// return 0;

	// long double p = 0.0l;
	// std::cout << "P(G | X, T ) = " << seq.pg_given_xtheta(g,x,theta) << std::endl;
	// std::cout << "px = " << seq.px_given_gtheta(x,g,theta) << std::endl;
	// std::cout << "p = " << exp(seq.px_given_gtheta(x,g,theta) + seq.pg(g)) << std::endl;
	// for(Genotype i : Genotype::enumerate_gts(2)){
	// 	long double px_i = seq.px_given_gtheta(x,i,theta);
	// 	std::cout << "G = " << i << "\tpx_i = " << px_i << std::endl;
	// 	p += exp(px_i + seq.pg(i));
	// }
	// std::cout << "divisor: " << p << std::endl;
	// std::cout << "result: " << exp(seq.px_given_gtheta(x,g,theta) + seq.pg(g)) / p << std::endl;
	// return 0;


	std::tuple<double> result = seq.start(.00001);
	std::cout << "Theta is: " << std::get<0>(result) << std::endl;
	return 0;
}