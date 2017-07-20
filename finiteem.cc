#include "finiteem.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <limits>

Finiteem::Finiteem(Pileupdata p, int ploidy) : plp(p), theta(std::make_tuple(0.1,std::map::map())),
	em(std::bind(&Finiteem::q_function, this, std::placeholders::_1), std::bind(&Finiteem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy){
	possible_gts = Genotype::enumerate_gts(ploidy);
}

Finiteem::Finiteem(Pileupdata p): Finiteem(p, 2){
}

Finiteem::Finiteem(std::string samfile, std::string refname, int ploidy) : plp(samfile, refname), theta(std::make_tuple(0.1,std::map::map())),
	em(std::bind(&Finiteem::q_function, this, std::placeholders::_1), std::bind(&Finiteem::m_function,this,std::placeholders::_1), theta),
	ploidy(ploidy){
	possible_gts = Genotype::enumerate_gts(ploidy);
}

Finiteem::Finiteem(std::string samfile, std::string refname) : Finiteem(samfile, refname, 2) {
}

theta_t Finiteem::start(double stop){
	return em.start(stop);
}

double Finiteem::q_function(theta_t theta){
	double likelihood = 0.0;
	for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
		likelihood += log(g.p_finite_alleles('A',0,std::get<0>(theta),std::get<1>(theta)));
	}
	return likelihood;
}

theta_t Finiteem::m_function(theta_t theta){
	pileupdata_t plpdata = plp.get_data();
	std::map<Genotype,int> gt_counts(possible_gts.begin(),possible_gts.end());
	std::map<char,std::map<Genotype,double>> t();
	std::vector<double> s(3,0.0);
	for (auto a: Genotype.alleles){
		t[*a] = gt_counts;
	}

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			const std::vector<char> &x = std::get<0>(*pos);
			const char &ref = std::get<3>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				std::vector<double> site_s = calc_s(x,*g,theta);
				double pg_x = pg_x_given_theta(*g,x,theta);
				for(size_t i = 0; i < s.size(); ++i){
					s[i] += pg_x * site_s[i];	
				}
				t[ref][*g] += pg_x;
			}
		}
	}

	// scale to prevent underflow
	// double smallest = smallest_nonzero(s);
	// if (smallest != 0){ //there SHOULD always be an s > 0.
	// 	std::transform(s.begin(),s.end(),s.begin(),[smallest](double d){ return d / smallest; });
	// }

	//quadratic formula
	double a = 3.0 * (s[0] + s[1] + s[2]);
	double b = - (3.0/2 * s[0] + s[1] + 5.0/2 * s[2]);
	double c = s[2] / 2;
	double epsilon_minus = (-b - sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
	// double epsilon_plus = (-b + sqrt(std::pow(b,2) - 4 * a * c))/(2 * a);
	double epsilon;
	if (epsilon_minus < 0){
		epsilon = 0;
	}
	else{
		epsilon = epsilon_minus;
	}
	return std::make_tuple(epsilon);
}

std::vector<double> Seqem::calc_s(std::vector<char> x, Genotype g, theta_t theta){ //TODO: make this generic
	std::vector<double> s(3,0.0);
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
double Seqem::pg_x_given_theta(const Genotype g,const std::vector<char> x,const theta_t theta){
	double px = px_given_gtheta(x,g,theta);
	return exp(px + pg(g));
}

//may be faster if we represent x as a map w/ char and counts, like gt?? we support this in plpdata
double Seqem::px_given_gtheta(const std::vector<char> x,const Genotype g,const theta_t theta){
	double px = 0.0;
	std::map<char, int> counts;

	for (std::vector<char>::const_iterator i = x.begin(); i != x.end(); ++i){
		counts[*i] += 1;
	}
	for (std::map<char,int>::iterator i = counts.begin(); i!= counts.end(); ++i){
		double pn = pn_given_gtheta(i->first,g,theta);
		if (pn == -std::numeric_limits<double>::infinity()){
			return -std::numeric_limits<double>::infinity();
		}
		else{
			px += i->second * pn;
			// px -= std::lgamma(i->second + 1);
		}
	}
	// px += std::lgamma(x.size()+1);
	return px;
}

//LOG SPACE
double Seqem::pn_given_gtheta(char n, Genotype g, theta_t theta){
	double epsilon = std::get<0>(theta);
	double p;
	// p = ((double)g.numbase(n))/g.getploidy()*(1.0-3.0*epsilon) + ((double)g.numnotbase(n))/g.getploidy()*epsilon;
	int numgt = g.numbase(n);
	if (numgt == 2){
		p = (1.0 - 3.0 * epsilon);
	}
	else if (numgt == 1){
		p = (0.5 - epsilon);
	}
	else{
		p = epsilon;
	}
	if (p == 0){
		return -std::numeric_limits<double>::infinity();
	}
	else if (p < 0){
		std::clog << "N=" << n << "\tGT=" << g << "\tPloidy=" << g.getploidy() << "\t#N=" << g.numbase(n) << "\t#!N=" << g.numnotbase(n) << "\tP1=" << ((double)g.numbase(n))/g.getploidy()*(1-3*epsilon) << "\tP2=" << ((double)g.numnotbase(n))/g.getploidy()*epsilon << "\tP=" << p << "\tTheta=" << theta << std::endl;
		throw std::runtime_error("p < 0 detected");
	}
	else{
		return log(p);
	}
}

double Seqem::pg(Genotype g){
	return log(1.0/possible_gts.size());
}

double Seqem::smallest_nonzero(std::vector<double> v){
	std::vector<double> sorted_v(v);
	std::sort(sorted_v.begin(),sorted_v.end());
	double smallest = 0;
	for (std::vector<double>::iterator i = sorted_v.begin(); i != sorted_v.end(); ++i){
		if(*i > smallest){
			return *i;
		}
	}
	return 0;
}

