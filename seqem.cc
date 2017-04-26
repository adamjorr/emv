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

long double Seqem::q_function(theta_t theta){
	long double p = 0;
	pileupdata_t plpdata = plp.get_data();

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(vector<pileuptuple_t>::iterator pos = *tid.begin(); pos != *tid.end(); ++pos){
			vector<char> x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				p += exp(pg_given_xtheta(*g, x, theta)) * (px_given_gtheta(x, *g, theta) + pg(*g));		
			}
		}
	}
	return p;
}

theta_t Seqem::m_function(theta_t theta){
	long double mu = std::get<0>(theta);
	pileupdata_t plpdata = plp.get_data();
	std::vector<long double> s(3,0); //TODO:make this generic, depends on ploidy

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(vector<pileuptuple_t>::iterator pos = *tid.begin(); pos != *tid.end(); ++pos){
			vector<char> x = std::get<0>(*pos);
			std::map<char,int> x_counts = std::get<1>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				std::vector<long double> site_s = calc_s(x,mu); //this would probably work better as a map: s[pn] += pg_given_xtheta + pn for each n
				for(int i = 0; i < s.length(); ++i){
					s[i] += pg_given_xtheta(x,mu) + site_s[i];	
				}
			}
		}
	}

	std::transform(s.begin(), s.end(), s.begin(), exp);

	long double mu_plus = (1.5 * s[0] + s[1] + 1.5 * s[2] + sqrt((1.5 * s[0] + s[1] - .5 * s[2])^2 - 2 * s[1] * s[2])) / (6 * s[0] + 6 * s[1] + 2 * s[2]);
	long double mu_minus = (1.5 * s[0] + s[1] + 1.5 * s[2] - sqrt((1.5 * s[0] + s[1] - .5 * s[2])^2 - 2 * s[1] * s[2])) / (6 * s[0] + 6 * s[1] + 2 * s[2]);
	return std::make_tuple(mu_minus > 0 ? std::min(mu_minus,mu_plus) : mu_plus);
}

static long double Seqem::calc_s(std::vector<char> x, Genotype g, theta_t theta){ //TODO: make this generic
	long double mu = std::get<0>(theta);
	std::vector<long double> s(3,0);
	for (std::std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		long double pn = pn_given_gtheta(*i, g, theta) 
		if (pn == log(1 - 3 * mu)){
			s[0]++;
		}
		else if (pn == log(.5 - mu)){
			s[1]++;
		}
		else if (pn == log(mu)){
			s[2]++;
		}
	}


}

long double Seqem::pg_given_xtheta(Genotype g, vector<char> x, theta_t theta){
	long double mu = std::get<0>(theta);
	long double p = px_given_gtheta(x,g,theta) + pg(g);
	for (std::vector<Genotype>::iterator i = possible_gts.begin(); i != possible_gts.end(); ++i){
		p -= px_given_gtheta(x,*i,theta) + pg(*i);
	}
	return p;
}

//may be faster if we represent x as a map w/ char and counts, like gt?? we support this in plpdata
long double Seqem::px_given_gtheta(vector<char> x, Genotype g, theta_t theta){
	long double p = 0;
	long double mu = std::get<0>(theta);
	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		p += pn_given_gtheta(*i, g, theta);
	}
	return p;
}

long double Seqem::pn_given_gtheta(char n, Genotype g, theta_t theta){
	return log(g.numbase(n)/g.getploidy()*(1-3*mu) + g.numnotbase(n)/g.getploidy()*mu))
}

long double Seqem::pg(Genotype g){
	return 1;
}
