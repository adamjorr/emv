#include "seqem.h"

Seqem::Seqem(Pileup plp, std::function<long double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta):
	em(q_function,m_function,theta), plp(plp) {}





