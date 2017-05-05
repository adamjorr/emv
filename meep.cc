#include "seqem.h"
#include "plpdata.h"
#include <string>
#include <iostream>

int main(int argc, char *argv[]){
	Seqem seq("testdata/test.sam","testdata/test.fa");
	std::cout << "Successfully obtained Seqem object";
	std::tuple<long double> result = seq.start(.001);
	std::cout << "Theta is: " << std::get<0l>(result) << std::endl;
}