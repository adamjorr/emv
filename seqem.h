#ifndef __MEEP_SEQEM_INCLUDED__
#define __MEEP_SEQEM_INCLUDED__

#include "pileup.h"
#include "em.h"

template<typename... T>
class Seqem{
protected:
	Pileup plp;
	EM em;
public:
	Seqem(std::string samfile, std::function<long double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta);
};




#endif
