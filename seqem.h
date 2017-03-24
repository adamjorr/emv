#ifndef __MEEP_SEQEM_INCLUDED__
#define __MEEP_SEQEM_INCLUDED__

#include "pileup.h"
#include "em.h"

class Seqem{
protected:
	Pileup plp;
	EM em;
public:
	long double q_function(std::vector<long double> theta);
	long double m_function(std::vector<long double> theta);
};




#endif
