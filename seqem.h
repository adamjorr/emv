#ifndef __MEEP_SEQEM_INCLUDED__
#define __MEEP_SEQEM_INCLUDED__

#include "pileup.h"
#include "plpdata.h"
#include "em.h"
#include "genotype.h"
#include <vector>

typedef std::tuple<double, double> theta_t; //mu parameter

class Seqem{
protected:
	Pileupdata plp;
	theta_t theta;
	EM<double> em;
	int ploidy;
	std::vector<Genotype> possible_gts;
public:
	Seqem(Pileupdata p, int ploidy);
	Seqem(Pileupdata p);
	Seqem(std::string samfile, std::string refname);
	Seqem(std::string samfile, std::string refname, int ploidy);
	theta_t start(double stop);
	double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
	static std::vector<double> calc_s(std::vector<char> x, Genotype g, theta_t theta);
	double pg_x_given_theta(Genotype g, std::vector<char> x, theta_t theta); //not log space
	double px_given_gtheta(std::vector<char> x, Genotype g, theta_t theta); // log space
	static double pn_given_gtheta(char n, Genotype g, theta_t theta); //log space
	double pg(Genotype g); //log space
	double smallest_nonzero(const std::vector<double> v);
};




#endif
