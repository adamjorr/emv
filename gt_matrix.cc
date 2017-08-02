#include "gt_matrix.h"
#include "popstatem.h"
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

GT_Matrix::GT_Matrix(std::string filename, int ploidy):
ploidy(ploidy), gts(Genotype::enumerate_gts(ploidy)), data(std::vector<double>(gts.size()),Genotype::alleles.size()){
	std::ifstream f;
	f.exceptions ( ifstream::failbit | ifstream::badbit );
	f.open(filename);
	std::string line;
	while(std::getline(f,line)){
		
	}
}

GT_Matrix::GT_Matrix(Pileupdata plpdata, int ploidy, theta_t theta, std::map<char,double> pi):
ploidy(ploidy), gts(Genotype::enumerate_gts(ploidy)), data(std::vector<double>(gts.size()),Genotype::alleles.size()){
	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			const std::vector<char> &x = std::get<0>(*pos);
			const char &ref = std::get<3>(*pos);
			for (std::vector<Genotype>::iterator g = gts.begin(); g != gts.end(); ++g){
				double pg_x = Popstatem::pg_x_given_theta(*g,x,theta,pi);
				double pdata = Popstatem::pdata_given_theta(x,theta,gts);
				this(ref,*g) += pg_x / pdata;
			}
		}
	}
}

std::vector<double>& GT_Matrix::operator[](size_t i){
	return data[i];
}

std::vector<double>& GT_Matrix::operator[](char allele){
	return data[std::distance(data.begin(),std::find(data.begin(), data.end(), allele))]
}

double& GT_Matrix::operator()(char allele, Genotype gt){
	row = this[allele];
	ptrdiff_t idx = std::distance(row.begin(),std::find(row.begin(), row.end(), gt));
	return row[idx];
}

double& GT_Matrix::operator()(size_t i, size_t j){
	row = this[i];
	return row[j];
}



