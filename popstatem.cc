#include "popstatem.h"
#include "meep_math.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <limits>

Popstatem::Popstatem(Pileupdata p, int ploidy) : plp(p), theta(std::make_tuple(0.1,std::map<char,double>(),1)),
	em(std::bind(&Popstatem::q_function, this, std::placeholders::_1), std::bind(&Popstatem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy){
	possible_gts = Genotype::enumerate_gts(ploidy);
	m = load_matrix();
}

Popstatem::Popstatem(Pileupdata p): Popstatem(p, 2){
}

Popstatem::Popstatem(std::string samfile, std::string refname, int ploidy) : plp(samfile, refname), theta(std::make_tuple(0.1,std::map<char, double>(),1)),
	em(std::bind(&Popstatem::q_function, this, std::placeholders::_1), std::bind(&Popstatem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy){
	possible_gts = Genotype::enumerate_gts(ploidy);
	m = load_matrix();
}

Popstatem::Popstatem(std::string samfile, std::string refname) : Popstatem(samfile, refname, 2) {
}

theta_t Popstatem::start(double stop){
	return em.start(stop);
}

double Popstatem::q_function(theta_t theta){
	double likelihood = 0.0;
	for (size_t i = 0; i < Genotype::alleles.size(); ++i){
		for (size_t j = 0; j < possible_gts.size(); ++j){
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){
				char allele = it->first;
				for (int k = 0; k < it->second; ++k){
					likelihood += m[i][j] * log(allele_alpha(allele,Genotype::alleles[i],std::get<2>(theta),std::get<0>(theta),std::get<1>(theta)) + k);
				}
			}
			for (int l = 0; l < ploidy; ++l){
				likelihood -= m[i][j] * log(ref_alpha(std::get<2>(theta), std::get<0>(theta) + l));
			}
		}
	}
	return likelihood;
}

theta_t Popstatem::m_function(theta_t theta){
	double th = meep_math::nr_root(std::bind(&Popstatem::dq_dtheta,this,std::placeholders::_1),std::bind(&Popstatem::ddq_dtheta,this,std::placeholders::_1),std::get<0>(theta));
	double w = meep_math::nr_root(std::bind(&Popstatem::dq_dw,this,std::placeholders::_1),std::bind(&Popstatem::ddq_dw,this,std::placeholders::_1),std::get<2>(theta));

	std::map<char,double> pi = std::get<1>(theta);
	double p = 0.0;
	for (size_t i = 0; i < Genotype::alleles.size()-1; ++i){
		char a = Genotype::alleles[i];
		auto qprime = std::bind(&Popstatem::dq_dpi,this,a,std::placeholders::_1);
		auto qprimeprime = std::bind(&Popstatem::ddq_dpi,this,a,std::placeholders::_1);
		double optimum = meep_math::nr_root(qprime,qprimeprime,pi[a]);
		pi[a] = optimum;
		p += optimum;
	}
	pi[Genotype::alleles.back()] = 1 - p;
	return std::make_tuple(th,pi,w);
}


//todo: figure out what to do with this function
GT_Matrix<10,4> Popstatem::load_matrix(){
	// for (pileupdata_t::iterator tid = plp.begin(); tid != plp.end(); ++tid){
	// 	for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
	// 		const std::vector<char> &x = std::get<0>(*pos);
	// 		const char &ref = std::get<3>(*pos);
	// 		// for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
	// 		// 	std::vector<double> site_s = calc_s(x,*g,theta);
	// 		// 	double pg_x = pg_x_given_theta(*g,x,theta);
	// 		// 	for(size_t i = 0; i < s.size(); ++i){
	// 		// 		s[i] += pg_x * site_s[i];	
	// 		// 	}
	// 		// 	t[ref][*g] += pg_x;
	// 		// }
	// 	}
	// }
	return GT_Matrix<10,4>();
}

double Popstatem::dq_dtheta(double th){
	std::map<char,double> pi = std::get<1>(theta);
	double refweight = std::get<2>(theta);
	double dq = 0.0;
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				for (int k = 0; k < it->second; ++k){ //add pi/(alpha + 0) twice for het; pi/(alpha + 1) for homozygote
					dq += m[i][j] * (pi[allele]/(allele_alpha(allele,Genotype::alleles[i],refweight,th,pi) + k));
				}
			}
			for (int l = 0; l < ploidy; ++l){ //subtract alpha, alpha + 1
				dq -= m[i][j] * (1.0 / ref_alpha(refweight, th) + l);	
			}
		}
	}
	return dq;
}

double Popstatem::ddq_dtheta(double th){
	double ddq = 0.0;
	std::map<char,double> pi = std::get<1>(theta);
	double refweight = std::get<2>(theta);
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (int l = 0; l < ploidy; ++l){ //add  1/alpha, 1/(alpha + 1)
				ddq += m[i][j] * (1.0 / pow(ref_alpha(refweight, th) + l,2));	
			}
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				for (int k = 0; k < it->second; ++k){ //add pi/(alpha + 0) twice for het; pi/(alpha + 1) for homozygote
					ddq -= m[i][j] * (pow(pi[allele],2)/pow(allele_alpha(allele,Genotype::alleles[i],refweight,th,pi) + k,2));
				}
			}
		}
	}
	return ddq;
}

double Popstatem::dq_dw(double w){
	double dq = 0.0;
	double th = std::get<0>(theta);
	std::map<char,double> pi = std::get<1>(theta);
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == Genotype::alleles[i]){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						dq += m[i][j] * (1.0/(allele_alpha(allele,(char)Genotype::alleles[i],w,th,pi) + k));
					}
				}
			}
		}
	}
	return dq;	
}

double Popstatem::ddq_dw(double w){
	double ddq = 0.0;
	double th = std::get<0>(theta);
	std::map<char,double> pi = std::get<1>(theta);
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == Genotype::alleles[i]){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						ddq -= m[i][j] * (1.0/pow(allele_alpha(allele,(char)Genotype::alleles[i],w,th,pi) + k,2));
					}
				}
			}
		}
	}
	return ddq;	
}

double Popstatem::dq_dpi(char a, double pi){
	double dq = 0.0;
	double th = std::get<0>(theta);
	double refweight = std::get<2>(theta);
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles - 1; ++i){ //for each reference base except one
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == a){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						dq += m[i][j] * ((th)/(allele_alpha(allele,(char)Genotype::alleles[i],refweight,th,pi) + k));
					}
				}
			}
		}
	}
	return dq;	
}

double Popstatem::ddq_dpi(char a, double pi){
	double ddq = 0.0;
	double th = std::get<0>(theta);
	double refweight = std::get<2>(theta);
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles - 1; ++i){ //for each reference base except one
		for (int j = 0; i < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == a){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						ddq -= m[i][j] * (pow(th,2)/pow(allele_alpha(allele,(char)Genotype::alleles[i],refweight,th,pi) + k,2));
					}
				}
			}
		}
	}
	return ddq;	
}

double Popstatem::allele_alpha(char allele, char ref, double ref_weight, double theta, std::map<char,double> pi){
	double w = (ref == allele ? ref_weight : 0);
	return theta * pi[allele] + w;
}

double Popstatem::allele_alpha(char allele, char ref, double ref_weight, double theta, double pi){
	double w = (ref == allele ? ref_weight : 0);
	return theta * pi + w;
}
double Popstatem::ref_alpha(double ref_weight, double theta){
	return ref_weight + theta;
}


