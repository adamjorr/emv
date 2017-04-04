#include "seqem.h"

Seqem::Seqem(std::string samfile, std::string refname, int ploidy) : em(q_function,m_function,theta), plp(samfile,refname), theta(), ploidy(ploidy) {
};

Seqem::Seqem(std::string samfile, std::string refname) : seqem(samfile, refname, 2) {
};

theta_t Seqem::start(long double stop){
	return em.start(stop);
}

long double q_function(theta_t theta){

}

theta_t m_function(theta_t theta){

}




