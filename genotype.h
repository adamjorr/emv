#ifndef __MEEP_GENOTYPE_INCLUDED__
#define __MEEP_GENOTYPE_INCLUDED__

class Genotype{
protected:
	map<char,int> gt;
	int ploidy;
public:
	Genotype(std::string gtstr);
	Genotype(map<char,int> gt);
	int numbase(char n);
	int numnotbase(char n);
	int getploidy();
};



#endif







