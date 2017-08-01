#ifndef __MEEP_GT_MATRIX_INCLUDED__
#define __MEEP_GT_MATRIX_INCLUDED__

#include <array>
#include "plpdata.h"
#include "popstatem.h"
#include "genotype.h"

template<size_t alleles, size_t gts>
class GT_Matrix{
protected:
	std::array<std::array<double,gts>,alleles> data;
public:
	double &operator[](size_t i);
	GT_Matrix(std::string filename);
	GT_Matrix(Pileupdata plpdata);
};

template<size_t alleles, size_t gts>
GT_Matrix<alleles,gts>::GT_Matrix(std::string filename){

}

template<size_t alleles, size_t gts>
GT_Matrix<alleles,gts>::GT_Matrix(Pileupdata plpdata, int ploidy, theta_t theta, std::map<char,double> pi){
	std::vector<Genotype,gts> possible_gts = Genotype::enumerate_gts(2);
	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			const std::vector<char> &x = std::get<0>(*pos);
			const char &ref = std::get<3>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				double pg_x = Seqem::pg_x_given_theta(*g,x,theta,pi);
			}
		}
	}
}

template<size_t alleles, size_t gts>
double &GT_Matrix<alleles,gts>::operator[](size_t i){
	return data[i];
}



#endif