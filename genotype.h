#ifndef __MEEP_GENOTYPE_INCLUDED__
#define __MEEP_GENOTYPE_INCLUDED__

class Genotype{
protected:
	map<char,int> gt;
	int ploidy;
	static std::vector<Genotype> enumerate_gts(std::vector<Genotype> v, int stopallele, int ploidy, std::string genotype);
	static const std::vector<char> alleles = {'A','T','G','C'};
public:
	Genotype(std::string gtstr);
	Genotype(map<char,int> gt);
	int numbase(char n);
	int numnotbase(char n);
	int getploidy();
	static std::vector<Genotype> enumerate_gts(int ploidy);
};



#endif







