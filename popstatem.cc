#include "popstatem.h"
#include "seqem.h"
#include "meep_math.h"
#include "gt_matrix.h"
#include "tuple_print.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <limits>
#include <boost/math/tools/roots.hpp>

typedef Popstatem::theta_t theta_t;

Popstatem::Popstatem(Pileupdata p, int ploidy) : plp(p), theta(std::make_tuple(0.1,Seqem::uniform_pi,1,0.1)),
	em(std::bind(&Popstatem::q_function, this, std::placeholders::_1), std::bind(&Popstatem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy), m(ploidy), possible_gts(Genotype::enumerate_gts(ploidy)){
	// possible_gts = Genotype::enumerate_gts(ploidy);
}

Popstatem::Popstatem(Pileupdata p): Popstatem(p, 2){
}

Popstatem::Popstatem(std::string samfile, std::string refname, int ploidy) : plp(samfile, refname), theta(std::make_tuple(0.1,Seqem::uniform_pi,1,0.1)),
	em(std::bind(&Popstatem::q_function, this, std::placeholders::_1), std::bind(&Popstatem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy), m(ploidy), possible_gts(Genotype::enumerate_gts(ploidy)){
	// possible_gts = Genotype::enumerate_gts(ploidy);
}

Popstatem::Popstatem(std::string samfile, std::string refname) : Popstatem(samfile, refname, 2) {
}

theta_t Popstatem::start(double stop){
	return em.start(stop);
}

double Popstatem::q_function(theta_t theta){
	double likelihood = 0.0;
	pileupdata_t plpdata = plp.get_data();

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			const std::vector<char> &x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				likelihood += pg_x_given_theta(*g,x,theta);
			}
		}
	}
	return likelihood;
}

theta_t Popstatem::m_function(theta_t theta){
	//optimize epsilon and update gt matrix m
	pileupdata_t plpdata = plp.get_data();
	std::vector<double> s(3,0.0); //TODO:make this generic, depends on ploidy
	GT_Matrix n;
	for (auto tid : plpdata){
		for(auto pos : tid){
			std::vector<char> x = std::get<0>(pos);
			char ref = std::get<3>(pos);
			double eps = std::get<3>(theta);
			std::map<char,double> pi = std::get<1>(theta);
			Seqem::increment_s(s, x, possible_gts, std::make_tuple(eps),pi);
			load_matrix(n,x,ref);
		}
	}
	double epsilon = Seqem::calc_epsilon(s);

	m = n; //this must be set before theta, w, and pi can be optimized
	std::cout << m << std::endl;

	//optimize theta, w, and pi
	double th = boost::math::tools::newton_raphson_iterate([this](const double& x){return std::make_tuple(dq_dtheta(x),ddq_dtheta(x));},std::get<0>(theta),0.0,10000.0,5);
	double w = boost::math::tools::newton_raphson_iterate([this](const double& x){return std::make_tuple(dq_dw(x),ddq_dw(x));},std::get<2>(theta),-1000.0,100000.0,5);

	std::map<char,double> pi = std::get<1>(theta);
	double p = 0.0;
	for (size_t i = 0; i < Genotype::alleles.size()-1; ++i){
		char a = Genotype::alleles[i];
		double optimum = boost::math::tools::newton_raphson_iterate([this,a](const double& x){return std::make_tuple(dq_dpi(a,x),ddq_dpi(a,x));},pi[a],0.0,1.0,5);
		std::cout << "allele:" << a << "\tpi:" << pi[a] <<"\toptimum:" << optimum  << std::endl;
		pi[a] = optimum;
		p += optimum;
	}
	pi[Genotype::alleles.back()] = 1 - p;

	return std::make_tuple(th,pi,w,epsilon);
}

void Popstatem::load_matrix(GT_Matrix &m, std::vector<char> x, char ref){
	for (auto g : possible_gts){
		double pg_x = pg_x_given_theta(g,x,theta);
		// double pdata = pdata_given_theta(x,theta,possible_gts);
		// m(ref,g) += (pg_x / pdata);
		m(ref,g) += pg_x;
	}
}

double Popstatem::dq_dtheta(double th){
	std::map<char,double> pi = std::get<1>(theta);
	double refweight = std::get<2>(theta);
	double dq = 0.0;
	int numalleles = Genotype::alleles.size();
	int numgts = possible_gts.size();
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; j < numgts; ++j){ //for each genotype
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
		for (int j = 0; j < numgts; ++j){ //for each genotype
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
		for (int j = 0; j < numgts; ++j){ //for each genotype
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
		for (int j = 0; j < numgts; ++j){ //for each genotype
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
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; j < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == a){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						dq += m[i][j] * th /(allele_alpha(allele,Genotype::alleles[i],refweight,th,pi) + k);
					}
				}
				// else{
				// 	for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
				// 		dq -= m[i][j] * th /(allele_alpha(allele,(char)Genotype::alleles[i],refweight,th,pi) + k) / 3;
				// 	}	
				// }
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
	for (int i = 0; i < numalleles; ++i){ //for each reference base
		for (int j = 0; j < numgts; ++j){ //for each genotype
			Genotype g = possible_gts[j];
			for (auto it = g.gt.begin(); it != g.gt.end(); ++it){ //for each base in genotype
				char allele = it -> first;
				if (allele == a){
					for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
						ddq -= m[i][j] * pow(th,2) / pow(allele_alpha(allele,Genotype::alleles[i],refweight,th,pi) + k,2);
					}
				}
				// else{
				// 	for (int k = 0; k < it->second; ++k){ //add 1/(alpha + 0) or +0 and +1 for homozygote
				// 		ddq -= m[i][j] * pow(th,2) / pow(allele_alpha(allele,(char)Genotype::alleles[i],refweight,th,pi) + k,2) / 9;
				// 	}	
				// }
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

double Popstatem::pg_x_given_theta(Genotype g, std::vector<char> x, theta_t theta){
	double e = std::get<3>(theta);
	std::map<char,double> p = std::get<1>(theta);
	return Seqem::pg_x_given_theta(g,x,std::make_tuple(e),p);
}

double Popstatem::pdata_given_theta(std::vector<char> x, theta_t theta, std::vector<Genotype> gts){
	double p = 0.0;
	for (auto g : gts){
		p += pg_x_given_theta(g, x, theta);
	}
	return p;
}



