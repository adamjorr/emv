#ifndef __MEEP_GT_MATRIX_INCLUDED__
#define __MEEP_GT_MATRIX_INCLUDED__

#include <array>
#include "plpdata.h"
#include "genotype.h"
#include "tuple_print.h"
#include <iostream>

//no reason for arrays over vectors:
//https://stackoverflow.com/questions/381621/using-arrays-or-stdvectors-in-c-whats-the-performance-gap
//cool faq on matrix classes:
//https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
//TODO: make a Matrix class, then a GT_Matrix class for allele+genotype based access, then make ctor a simple function in popstatem
class GT_Matrix{
protected:
	int ploidy;
	std::vector<Genotype> gts;
	std::vector<std::vector<double>> data;
public:
	friend std::ostream& operator<<(std::ostream& os, const GT_Matrix);
	std::vector<double> operator[](size_t i) const {return data[i];};
	std::vector<double> operator[](int i) const {return data[i];};
	std::vector<double> &operator[](size_t i) {return data[i];};
	std::vector<double> &operator[](int i) {return data[i];};
	const std::vector<double> &operator[](char allele) const;
	double operator()(char allele, Genotype gt) const;
	double operator()(size_t i, size_t j) const {return data[i][j];};
	double &operator()(char allele, Genotype gt);
	double &operator()(size_t i, size_t j) {return data[i][j];};
	GT_Matrix();
	GT_Matrix(int ploidy);
	GT_Matrix(std::string filename, int ploidy);
	// GT_Matrix(Pileupdata plpdata, int ploidy, theta_t theta, std::map<char,double> pi);
};

std::ostream& operator<<(std::ostream& os, const GT_Matrix);

#endif