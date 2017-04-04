#ifndef __MEEP_SEQEM_INCLUDED__
#define __MEEP_SEQEM_INCLUDED__

#include "pileup.h"
#include "em.h"

typedef std::tuple<long double> theta_t;

class Seqem{
protected:
	Pileupdata plp;
	EM em;
	theta_t theta;
	int ploidy;
public:
	Seqem(std::string samfile, std::string refname);
	Seqem(std::string samfile, std::string refname, int ploidy);
	theta_t start(long double stop);
	long double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
	long double pgt_given_xtheta(vector<char> x, theta_t theta);
	long double px_given_ztheta(vector<char> gt, )
};




#endif
