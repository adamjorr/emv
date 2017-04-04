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
public:
	Seqem(std::string samfile, std::string refname);
	theta_t start(long double stop);
	long double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
};




#endif
