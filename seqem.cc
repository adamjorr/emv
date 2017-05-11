#include "seqem.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>

Seqem::Seqem(std::string samfile, std::string refname, int ploidy) : 
	plp(samfile,refname), theta(std::make_tuple(0)),
	em(std::bind(&Seqem::q_function, this, std::placeholders::_1), std::bind(&Seqem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy) {
	possible_gts = Genotype::enumerate_gts(ploidy);
};

Seqem::Seqem(std::string samfile, std::string refname) : Seqem(samfile, refname, 2) {
};

theta_t Seqem::start(long double stop){
	return em.start(stop);
}

long double Seqem::q_function(theta_t theta){
	long double p = 0;
	pileupdata_t plpdata = plp.get_data();

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			std::vector<char> x = std::get<0>(*pos);
			// std::cout << "X=" << x << std::endl;
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				// p += exp(pg_given_xtheta(*g, x, theta)) * (px_given_gtheta(x, *g, theta) + pg(*g));
				long double a = pg_given_xtheta(*g, x, theta);
				if (a != 0){
					long double b = px_given_gtheta(x, *g, theta);
					long double c = pg(*g);
					p += a * (b + c);
					// std::cout << "P(G | X, T) = " << a << "\tG=" << *g << "\tT=" << theta << std::endl;
					// std::cout << "P(X | G, T) = " << b << std::endl;
					// std::cout << "P(G) = " << c << "\tG=" << *g << "\tT=" << theta << std::endl;
				}
				// if (errno){
				// 	std::cout << "P(G | X, T) = " << a << std::endl;
				// 	std::cout << "P(X | G, T) = " << b << std::endl;
				// 	std::cout << "P(G) = " << c << std::endl;
				// 	std::cout << "p = " << p << std::endl;
				// 	// std::cout << std::strerror(errno) << std::endl;
				// 	throw std::runtime_error(std::strerror(errno));
				// }
			}
		}
	}
	// std::cout << "Likelihood returned is: " << p << std::endl;
	return p;
}

theta_t Seqem::m_function(theta_t theta){
	pileupdata_t plpdata = plp.get_data();
	std::vector<long double> s(3,0); //TODO:make this generic, depends on ploidy

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			std::vector<char> x = std::get<0>(*pos);
			std::map<char,int> x_counts = std::get<1>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				std::vector<long double> site_s = calc_s(x,*g,theta); //this would probably work better as a map: s[pn] += pg_given_xtheta + pn for each n
				long double pg = pg_given_xtheta(*g,x,theta);
				if (pg != 0){
					for(int i = 0; i < s.size(); ++i){
						s[i] += pg * site_s[i];	
					}
				}
			}
		}
	}

	// scale to prevent underflow
	long double smallest = smallest_nonzero(s);
	if (smallest == 0){ //there SHOULD always be an s > 0.
		std::clog << "s: " << s << std::endl;
		throw std::logic_error("All s == 0, this should never happen.");
	}
	// std::clog << "s: " << s;
	std::transform(s.begin(),s.end(),s.begin(),[smallest](long double d){ return d / smallest; });
	// std::clog << " => S: " << s << std::endl;

	long double mu_plus = (1.5 * s[0] + s[1] + 1.5 * s[2] + sqrt(std::pow(1.5 * s[0] + s[1] - .5 * s[2],2) - 2 * s[1] * s[2])) / (6 * s[0] + 6 * s[1] + 2 * s[2]);
	long double mu_minus = (1.5 * s[0] + s[1] + 1.5 * s[2] - sqrt(std::pow(1.5 * s[0] + s[1] - .5 * s[2],2) - 2 * s[1] * s[2])) / (6 * s[0] + 6 * s[1] + 2 * s[2]);
	// if(errno){
	// 	std::clog << "Error calculating mu+ and mu-" << std::endl;
	// 	throw std::runtime_error(std::strerror(errno));
	// }
	// std::clog << "mu+: " << mu_plus << "\t mu-: " << mu_minus << std::endl;
	return std::make_tuple(mu_minus > 0 ? std::min(mu_minus,mu_plus) : mu_plus);
}

std::vector<long double> Seqem::calc_s(std::vector<char> x, Genotype g, theta_t theta){ //TODO: make this generic
	long double mu = std::get<0>(theta);
	std::vector<long double> s(3,0);
	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		long double pn = pn_given_gtheta(*i, g, theta);
		if (pn == (1 - 3 * mu)){
			s[0]++;
		}
		else if (pn == (.5 - mu)){
			s[1]++;
		}
		else if (pn == (mu)){
			s[2]++;
		}
	}
	if(errno){
		std::clog << "Error calculating s" << std::endl;
		throw std::runtime_error(std::strerror(errno));
	}
	return s;
}

//RESULT NOT IN LOG SPACE
long double Seqem::pg_given_xtheta(Genotype g, std::vector<char> x, theta_t theta){
	long double px = px_given_gtheta(x,g,theta);
	if (px==0){
		return px;
	}
	long double p = exp(px + pg(g));
	long double division_p = 0;
	// std::clog << "P(X | G, theta) = " << px_given_gtheta(x,g,theta) << std::endl; 
	// std::clog << "P(X | G, theta) * P(G) = " << p << std::endl;
	for (std::vector<Genotype>::iterator i = possible_gts.begin(); i != possible_gts.end(); ++i){
		long double px_i = px_given_gtheta(x,*i,theta);
		if (px_i != 0){
			division_p += exp(px_i + pg(*i));
		}
		// if(errno){
		// 	std::clog << "Error calculating pg_given_xtheta" << std::endl;
		// 	std::clog << "p: " << p << std::endl;
		// 	throw std::runtime_error(std::strerror(errno));
		// }
	}
	return p / division_p;
}

//may be faster if we represent x as a map w/ char and counts, like gt?? we support this in plpdata
long double Seqem::px_given_gtheta(std::vector<char> x, Genotype g, theta_t theta){
	long double px = 0;
	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		long double pn = pn_given_gtheta(*i, g, theta);
		if (pn > 0){
			px += log(pn);
		}
		// if(errno){
		// 	std::clog << "Error calculating px_given_gtheta" << std::endl;
		// 	std::clog << "P(N): " << pn << "\tP(X|G,theta): " << px << std::endl;
		// 	throw std::runtime_error(std::strerror(errno));
		// }
	}
	// if (px == 0){
	// 	std::clog << "Error calculating px_given_gtheta" << std::endl;
	// 	std::clog << "X=" << x << "\tG=" << g << "\tP(X|G,theta): " << px << std::endl;
	// 	throw std::runtime_error("px_given_gtheta should never be 0.");
	// }
	return px;
}

//THIS IS NOT IN LOG SPACE
long double Seqem::pn_given_gtheta(char n, Genotype g, theta_t theta){
	long double mu = std::get<0>(theta);
	long double p = ((long double)g.numbase(n))/g.getploidy()*(1-3*mu) + ((long double)g.numnotbase(n))/g.getploidy()*mu;
	// std::clog << "N=" << n << "\tGT=" << g << "\tPloidy=" << g.getploidy() << "\t#N=" << g.numbase(n) << "\t#!N=" << g.numnotbase(n) << "\tP1=" << ((long double)g.numbase(n))/g.getploidy()*(1-3*mu) << "\tP2=" << ((long double)g.numnotbase(n))/g.getploidy()*mu << "\tP=" << p << "\tTheta=" << theta << std::endl;
	// p = (p <= 0 ? 0 : log(p));
	// if(errno){
	// 	std::clog << "Error calculating pn_given_gtheta" << std::endl;
	// 	std::clog << "N=" << n << "\tGT=" << g << "\t#N=" << g.numbase(n) << "\t#!N=" << g.numnotbase(n) << "\tP=" << p << "\tTheta=" << theta << std::endl;
	// 	throw std::runtime_error(std::strerror(errno));
	// }
	return p;
}

long double Seqem::pg(Genotype g){
	return log(1.0l/possible_gts.size());
}

long double Seqem::smallest_nonzero(std::vector<long double> v){
	std::vector<long double> sorted_v(v);
	std::sort(sorted_v.begin(),sorted_v.end());
	long double smallest = 0;
	for (std::vector<long double>::iterator i = sorted_v.begin(); i != sorted_v.end(); ++i){
		if(*i > smallest){
			return *i;
		}
	}
	return 0;
}

