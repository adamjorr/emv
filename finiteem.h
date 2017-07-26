#ifndef __MEEP_FINITEEM_INCLUDED__
#define __MEEP_FINITEEM_INCLUDED__

#include "plpdata.h"
#include "em.h"
#include "genotype.h"
#include <vector>

typedef std::tuple<double, std::map<char, int>, double> theta_t; //theta, pi, refweight

class Seqem{
protected:
	Pileupdata plp;
	theta_t theta;
	EM<double> em;
	int ploidy;
	std::vector<Genotype> possible_gts;
	GT_Matrix<ploidy> m;
public:
	Finiteem(Pileupdata p, int ploidy);
	Finiteem(Pileupdata p);
	Finiteem(std::string samfile, std::string refname);
	Finiteem(std::string samfile, std::string refname, int ploidy);
	theta_t start(double stop);
	double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
	static double allele_alpha(char allele, char ref, double ref_weight, double theta, std::map<char,double> pi);
	static double ref_alpha(double ref_weight, double theta);
};

template<int ploidy>
using GT_Matrix = std::array<std::array<int,Genotype.enumerate_gts(ploidy).size()>,Genotype.alleles.size()>;

#endif
