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
	std::vector<double> s(3,0.0); //TODO:make this generic, depends on ploidy

	for (pileupdata_t::iterator tid = plpdata.begin(); tid != plpdata.end(); ++tid){
		for(std::vector<pileuptuple_t>::iterator pos = tid->begin(); pos != tid->end(); ++pos){
			const std::vector<char> &x = std::get<0>(*pos);
			for (std::vector<Genotype>::iterator g = possible_gts.begin(); g != possible_gts.end(); ++g){
				std::vector<double> site_s = calc_s(x,*g,theta);
				double pg_x = pg_x_given_theta(*g,x,theta);
				for(size_t i = 0; i < s.size(); ++i){
					s[i] += pg_x * site_s[i];	
				}
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