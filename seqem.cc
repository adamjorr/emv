#include "seqem.h"

Seqem::Seqem(std::string samfile, std::string refname) : em(q_function,m_function,theta), plp(samfile,refname), theta() {
};

theta_t Seqem::start(long double stop){
	return em.start(stop);
}

long double q_function(theta_t theta){

}

theta_t m_function(theta_t theta){

}




