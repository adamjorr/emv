#include "seqem.h"
#include "plpdata.h"
#include <string>
#include <iostream>
#include <cmath>

int main(int argc, char *argv[]){
	// Seqem seq("testdata/test.sam","testdata/test.fa");
	std::clog.precision(20);
	Seqem seq("foo.sam","testdata/test.fa");

	// Genotype g("CC");
	// std::vector<char> x(26,'C');
	// std::tuple<long double> theta = std::make_tuple(0.0l);
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


	std::tuple<long double> result = seq.start(.00001);
	std::cout << "Theta is: " << std::get<0>(result) << std::endl;
	return 0;
}