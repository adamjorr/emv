#ifndef __MEEP_GT_MATRIX_INCLUDED__
#define __MEEP_GT_MATRIX_INCLUDED__

#include <array>
#include "plpdata.h"
#include "popstatem.h"
#include "genotype.h"

//no reason for arrays over vectors:
//https://stackoverflow.com/questions/381621/using-arrays-or-stdvectors-in-c-whats-the-performance-gap
//cool faq on matrix classes:
//https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
class GT_Matrix{
protected:
	std::vector<Genotype> gts;
	std::vector<std::vector<double>> data;
	int ploidy;
public:
	std::vector<double> &operator[](size_t i);
	std::vector<double> &operator[](char allele);
	double &operator()(char allele, Genotype gt);
	double &operator()(size_t i, size_t j);
	GT_Matrix(std::string filename);
	GT_Matrix(Pileupdata plpdata);
};

#endif