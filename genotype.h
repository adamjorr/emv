#ifndef __MEEP_GENOTYPE_INCLUDED__
#define __MEEP_GENOTYPE_INCLUDED__

#include <map>
#include <vector>
#include <string>
#include <iostream>

class Genotype{
protected:
	int ploidy;
	static void enumerate_gts(std::vector<Genotype> &v, int stopallele, unsigned int ploidy, std::string genotype);
public:
	std::map<char,int> gt;
	Genotype(std::string gtstr);
	Genotype(std::map<char,int> gt);
	double p_finite_alleles(char ref, double ref_weight, double theta, std::map<char,double> pi);
	int numbase(char n);
	int numnotbase(char n);
	int getploidy();
	std::string to_string() const;
	static std::vector<Genotype> enumerate_gts(int ploidy);
	static const std::vector<char> alleles;
};

std::ostream& operator<<(std::ostream& os, const Genotype);

#endif







