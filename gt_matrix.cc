#include "gt_matrix.h"
#include <algorithm>
#include <iterator>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

GT_Matrix::GT_Matrix(): GT_Matrix(2) {}

GT_Matrix::GT_Matrix(int ploidy): ploidy(ploidy), gts(Genotype::enumerate_gts(ploidy)), data(Genotype::alleles.size(),std::vector<double>(gts.size(),0.0)) {
}

GT_Matrix::GT_Matrix(std::string filename, int ploidy):
ploidy(ploidy), gts(Genotype::enumerate_gts(ploidy)), data(Genotype::alleles.size(),std::vector<double>(gts.size())){
	std::ifstream f;
	f.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	f.open(filename);
	int i = 0;
	for(std::string line; std::getline(f,line);){
		std::istringstream row(line);
		int j = 0;
		for(std::string val; std::getline(row,val,'\t');){
			data[i][j] = std::stod(val);
		}
	}
	f.close();
}

// GT_Matrix::GT_Matrix(Pileupdata plpdata, int ploidy, theta_t theta):
// ploidy(ploidy), gts(Genotype::enumerate_gts(ploidy)), data(std::vector<double>(gts.size()),Genotype::alleles.size()){
// 	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
// 		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
// 			const std::vector<char> &x = std::get<0>(*pos);
// 			const char &ref = std::get<3>(*pos);
// 			for (std::vector<Genotype>::iterator g = gts.begin(); g != gts.end(); ++g){
// 				double pg_x = Popstatem::pg_x_given_theta(*g,x,theta);
// 				double pdata = Popstatem::pdata_given_theta(x,theta,gts);
// 				this(ref,*g) += pg_x / pdata;
// 			}
// 		}
// 	}
// }

const std::vector<double>& GT_Matrix::operator[](char allele) const{
	return data[std::distance(Genotype::alleles.begin(),std::find(Genotype::alleles.begin(), Genotype::alleles.end(), allele))];
}

double GT_Matrix::operator()(char allele, Genotype gt) const{
	std::vector<double> row = data[std::distance(Genotype::alleles.begin(),std::find(Genotype::alleles.begin(), Genotype::alleles.end(), allele))];
	ptrdiff_t idx = std::distance(gts.begin(),std::find(gts.begin(), gts.end(), gt));
	return row[idx];
}

double& GT_Matrix::operator()(char allele, Genotype gt){
	std::vector<double>& row = data[std::distance(Genotype::alleles.begin(),std::find(Genotype::alleles.begin(), Genotype::alleles.end(), allele))];
	ptrdiff_t idx = std::distance(gts.begin(),std::find(gts.begin(), gts.end(), gt));
	return row[idx];
}

std::ostream& operator<<(std::ostream& os, const GT_Matrix m){
	os << "- matrix -" << std::endl;
	for(int i = 0; i < Genotype::alleles.size(); ++i){
		os << m[i] << std::endl;
	}
	return os << "----------" << std::endl;
}

