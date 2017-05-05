#include "genotype.h"

Genotype::Genotype(std::string gtstr) : gt(), ploidy(0) {
	for(auto it = gtstr.begin(); it != gtstr.end(); ++it){
		++gt[*it];
		ploidy+=1;
	}
}

Genotype::Genotype(std::map<char,int> gt) : gt(gt), ploidy(0) {
	for(auto it = gt.begin(); it!= gt.end(); ++it){
		ploidy += it->second;
	}
}

int Genotype::numbase(char n){
	return gt[n];
}

int Genotype::numnotbase(char n){
	return ploidy - gt[n];
}

int Genotype::getploidy(){
	return ploidy;
}

//check out http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
void Genotype::enumerate_gts(std::vector<Genotype> v, int stopallele, unsigned int ploidy, std::string genotype){
	if (genotype.length()==ploidy){
		v.push_back(Genotype(genotype));
	}
	else{
		for (int i = 0; i < stopallele; ++i){
			std::string s(&alleles[i]);
			s.append(genotype);
			enumerate_gts(v,i+1, ploidy, s);
		}
	}
}

std::vector<Genotype> Genotype::enumerate_gts(int ploidy){
	std::vector<Genotype> v;
	enumerate_gts(v,alleles.size(),ploidy,"");
	return v;
}

const std::vector<char> Genotype::alleles = {'A','T','G','C'};




