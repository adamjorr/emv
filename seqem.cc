#include "seqem.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <limits>

Seqem::Seqem(std::string samfile, std::string refname, int ploidy) : 
	plp(samfile,refname), theta(std::make_tuple(0.1l)),
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
	long double p = 0.0l;
	long double likelihood = 0.0l;
	pileupdata_t plpdata = plp.get_data();

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			std::vector<char> x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				long double a = pg_given_xtheta(*g, x, theta);
				long double b = px_given_gtheta(x, *g, theta);
				long double c = pg(*g);
				p += a * (b + c);
				likelihood += b + c;
			}
		}
	}
	std::clog << likelihood << std::endl;
	return p;
}

theta_t Seqem::m_function(theta_t theta){
	pileupdata_t plpdata = plp.get_data();
	std::vector<long double> s(3,0.0l); //TODO:make this generic, depends on ploidy

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			std::vector<char> x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				std::vector<long double> site_s = calc_s(x,*g,theta);
				long double pg = pg_given_xtheta(*g,x,theta);
				for(int i = 0; i < s.size(); ++i){
					s[i] += pg * site_s[i];	
				}
			}
		}
	}

	// scale to prevent underflow
	long double smallest = smallest_nonzero(s);
	if (smallest != 0){ //there SHOULD always be an s > 0.
		std::transform(s.begin(),s.end(),s.begin(),[smallest](long double d){ return d / smallest; });
	}

	//quadratic formula
	long double a = (3l * s[0] + 3l * s[1] + s[2]);
	long double b = - (3l/2 * s[0] + s[1] + 3l/2 * s[2]);
	long double c = s[2] / 2;
	long double mu_minus = (-b - sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
	long double mu_plus = (-b + sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);

	//substitute
	// long double mu_minus = ((3l/2 * s[0] + s[1] + 3l/2 * s[2]) - sqrt(std::pow(- (3l/2 * s[0] + s[1] + 3l/2 * s[2]),2) - 4 * (3l * s[0] + 3l * s[1] + s[2]) * s[2] / 2))/(2 * (3l * s[0] + 3l * s[1] + s[2]));
	// long double mu_plus = ((3l/2 * s[0] + s[1] + 3l/2 * s[2]) + sqrt(std::pow(- (3l/2 * s[0] + s[1] + 3l/2 * s[2]),2) - 4 * (3l * s[0] + 3l * s[1] + s[2]) * s[2] / 2))/(2 * (3l * s[0] + 3l * s[1] + s[2]));
	long double mu = (mu_minus > 0 ? std::min(mu_minus,mu_plus) : mu_plus);
	return std::make_tuple(mu);
}

std::vector<long double> Seqem::calc_s(std::vector<char> x, Genotype g, theta_t theta){ //TODO: make this generic
	long double mu = std::get<0>(theta);
	std::vector<long double> s(3,0.0l);
	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		int numgt = g.numbase(*i);
		if (numgt == 2){
			s[0]++;
		}
		else if (numgt == 1){
			s[1]++;
		}
		else if (numgt == 0){
			s[2]++;
		}
	}
	return s;
}

//RESULT NOT IN LOG SPACE
long double Seqem::pg_given_xtheta(Genotype g, std::vector<char> x, theta_t theta){
	long double px = px_given_gtheta(x,g,theta);
	long double p = exp(px + pg(g));
	if (p == 0){
		return 0;
	}
	long double division_p = 0.0l;
	for (std::vector<Genotype>::iterator i = possible_gts.begin(); i != possible_gts.end(); ++i){
		long double px_i = px_given_gtheta(x,*i,theta);
		division_p += exp(px_i + pg(*i));
	}
	if (division_p == 0){
		std::clog << "Error calculating pg_given_xtheta" << std::endl;
		std::clog << "G=" << g << std::endl;
		std::clog << "X=" << x << std::endl;
		std::clog << "T=" << theta << std::endl;
		std::clog << "p: " << p << std::endl;
		throw std::runtime_error("pg_given_xtheta");
	}
	return p / division_p;
}

//may be faster if we represent x as a map w/ char and counts, like gt?? we support this in plpdata
long double Seqem::px_given_gtheta(std::vector<char> x, Genotype g, theta_t theta){
	long double px = 0.0l;
	std::map<char, int> counts;

	std::cout << "X = " << x << std::endl;
	std::cout << "G = " << g << std::endl;
	std::cout << "T = " << theta << std::endl;

	for (std::vector<char>::iterator i = x.begin(); i != x.end(); ++i){
		counts[*i] += 1;
	}

	std::cout << "Counts = " << counts << std::endl;

	for (std::map<char,int>::iterator i = counts.begin(); i!= counts.end(); ++i){
		long double pn = pn_given_gtheta(i->first,g,theta);
		if (pn == -std::numeric_limits<long double>::infinity()){
			return -std::numeric_limits<long double>::infinity();
		}
		else{
			px += i->second * pn_given_gtheta(i->first, g, theta);
			px -= std::lgamma(i->second + 1);
		}
	}
	px += std::lgamma(x.size()+1);

	std::cout << "P(X | G, T) = " << px << std::endl;

	return px;
}

//LOG SPACE
long double Seqem::pn_given_gtheta(char n, Genotype g, theta_t theta){
	long double mu = std::get<0>(theta);
	long double p;
	// p = ((long double)g.numbase(n))/g.getploidy()*(1.0l-3.0l*mu) + ((long double)g.numnotbase(n))/g.getploidy()*mu;
	int numgt = g.numbase(n);
	if (numgt == 2){
		p = (1.0l - 3.0l * mu);
	}
	else if (numgt == 1){
		p = (0.5l - mu);
	}
	else{
		p = mu;
	}
	if (p == 0){
		return -std::numeric_limits<long double>::infinity();
	}
	else if (p < 0){
		std::clog << "N=" << n << "\tGT=" << g << "\tPloidy=" << g.getploidy() << "\t#N=" << g.numbase(n) << "\t#!N=" << g.numnotbase(n) << "\tP1=" << ((long double)g.numbase(n))/g.getploidy()*(1-3*mu) << "\tP2=" << ((long double)g.numnotbase(n))/g.getploidy()*mu << "\tP=" << p << "\tTheta=" << theta << std::endl;
		throw std::runtime_error("p < 0 detected");
	}
	else{
		return log(p);
	}
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

