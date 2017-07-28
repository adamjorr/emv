#ifndef __MEEP_POPSTATEM_INCLUDED__
#define __MEEP_POPSTATEM_INCLUDED__

#include "plpdata.h"
#include "em.h"
#include "genotype.h"
#include <vector>

typedef std::tuple<double, std::map<char, double>, double> theta_t; //theta, pi, refweight

template<int alleles, int gts>
using GT_Matrix = std::array<std::array<double,gts>,alleles>;
// using GT_Matrix = std::array<std::array<int,Genotype::enumerate_gts(ploidy).size()>,Genotype::alleles.size()>;

class Popstatem{
protected:
	Pileupdata plp;
	theta_t theta;
	EM<double,std::map<char,double>,double> em;
	int ploidy;
	std::vector<Genotype> possible_gts;
	GT_Matrix<10,4> m;
public:
	Popstatem(Pileupdata p, int ploidy);
	Popstatem(Pileupdata p);
	Popstatem(std::string samfile, std::string refname);
	Popstatem(std::string samfile, std::string refname, int ploidy);
	theta_t start(double stop);
	double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
	GT_Matrix<10,4> load_matrix();
	double dq_dtheta(double th);
	double ddq_dtheta(double th);
	double dq_dw(double w);
	double ddq_dw(double w);
	double dq_dpi(char a, double pi);
	double ddq_dpi(char a, double pi);
	static double allele_alpha(char allele, char ref, double ref_weight, double theta, std::map<char,double> pi);
	static double allele_alpha(char allele, char ref, double ref_weight, double theta, double pi);
	static double ref_alpha(double ref_weight, double theta);
};

#endif
