#include "seqem.h"
#include "plpdata.h"
#include <string>
#include <iostream>
#include <cmath>

int main(int argc, char *argv[]){
	std::clog.precision(15);

	Seqem seq("foo.sam","testdata/test.fa");


	// std::vector<char> x(9,'C');
	// x.push_back('A');
	// Pileupdata data(x);
	// Seqem seq(data);


	std::tuple<double> result = seq.start(.00001);
	std::cout << "Theta is: " << std::get<0>(result) << std::endl;
	return 0;
}