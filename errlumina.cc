#include "samio.h"
#include "mismatchfinder.h"
#include <string>
#include <iostream>

int main(int argc, char *argv[]){
	SamReader reader("testdata/test.sam");
	MismatchFinder mismatcher(&reader,"testdata/test.fa");
	mismatcher.dump_locations();
}