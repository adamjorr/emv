#include "genotype.h"

Genotype::Genotype(std::string gtstr) : gt(), ploidy(0) {
	for(auto it = gtstr.begin(); it != gtstr.end(); ++it){
		++gt[*it];
		ploidy+=1;
	}
}

Genotype::Genotype(map<char,int> gt) : gt(gt), ploidy(0) {
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





