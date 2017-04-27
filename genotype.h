#ifndef __MEEP_GENOTYPE_INCLUDED__
#define __MEEP_GENOTYPE_INCLUDED__

#include <map>
#include <vector>
#include <string>

class Genotype{
protected:
	std::map<char,int> gt;
	int ploidy;
	static void enumerate_gts(std::vector<Genotype> v, int stopallele, int ploidy, std::string genotype);
	static const std::vector<char> alleles;
public:
	Genotype(std::string gtstr);
	Genotype(std::map<char,int> gt);
	int numbase(char n);
	int numnotbase(char n);
	int getploidy();
	static std::vector<Genotype> enumerate_gts(int ploidy);
};



#endif







