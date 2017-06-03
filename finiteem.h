#ifndef __MEEP_FINITEEM_INCLUDED__
#define __MEEP_FINITEEM_INCLUDED__

#include "plpdata.h"
#include "em.h"
#include "genotype.h"
#include <vector>

typedef std::tuple<double, std::map<char, int> > theta_t; //theta, pi

class Seqem{
protected:
	Pileupdata plp;
	theta_t theta;
	EM<double> em;
	int ploidy;
	std::vector<Genotype> possible_gts;
public:
	Finiteem(Pileupdata p, int ploidy);
	Finiteem(Pileupdata p);
	Finiteem(std::string samfile, std::string refname);
	Finiteem(std::string samfile, std::string refname, int ploidy);
	theta_t start(double stop);
	double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
};




#endif
