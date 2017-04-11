#include "seqem.h"
#include <cmath>

Seqem::Seqem(std::string samfile, std::string refname, int ploidy) : em(q_function,m_function,theta), plp(samfile,refname), theta(), ploidy(ploidy) {
	possible_gts = Genotype.enumerate_gts(ploidy);
};

Seqem::Seqem(std::string samfile, std::string refname) : seqem(samfile, refname, 2) {
};

theta_t Seqem::start(long double stop){
	return em.start(stop);
}

long double q_function(theta_t theta){
	long double p = 0;
	pileupdata_t plpdata = plp.get_data();

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(vector<pileuptuple_t>::iterator pos = *tid.begin(); pos != *tid.end(); ++pos){
			vector<char> x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				p += pg_given_xtheta(*g, x, theta) + log(px_given_gtheta(x, *g, theta) + pg(*g));		
			}
		}
	}
	return p;
}

theta_t m_function(theta_t theta){
	
}

long double pg_given_xtheta(Genotype g, vector<char> x, theta_t theta){
	long double mu = std::get<0>(theta);
	long double p = px_given_gtheta(x,g,theta) + pg(g);
	for (std::vector<Genotype>::iterator i = possible_gts.begin(); i != possible_gts.end(); ++i){
		p -= px_given_gtheta(x,*i,theta) + pg(*i);
	}
	return p;
}

//may be faster if we represent x as a map w/ char and counts, like gt?? we support this in plpdata
long double px_given_gtheta(vector<char> x, Genotype g, theta_t theta){
	long double p = 0;
	long double mu = std::get<0>(theta);
	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		p += log(g.numbase(*i)/g.getploidy()*(1-3*mu) + g.numnotbase(*i)/g.getploidy()*mu);
	}
	return p;
}

long double pg(Genotype g){

}
