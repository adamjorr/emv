#include "popstatem.h"
#include "plpdata.h"
#include <string>
#include <iostream>
#include <cmath>
#include "meep_math.h"

int main(int argc, char *argv[]){
	std::clog.precision(15);

	// Seqem seq("foo.sam","testdata/test.fa");
	// std::vector<char> x(9,'C');
	// x.push_back('A');
	// Pileupdata data(x);
	// Seqem seq(data);

	Popstatem seq("foo.sam","testdata/test.fa");
	std::tuple<double, std::map<char,double>, double, double> result = seq.start(.00001);
	std::cout << "Theta is: " << result << std::endl;

	// std::cout << meep_math::nr_root([](double x){return -4 * x + 3;}, [](double x){return -4;}, 10 ) << std::endl;


	return 0;
}