#ifndef __MEEP_POPSTATEM_INCLUDED__
#define __MEEP_POPSTATEM_INCLUDED__

#include "plpdata.h"
#include "em.h"
#include "genotype.h"
#include "gt_matrix.h"
#include <vector>

// template<int alleles, int gts>
// using GT_Matrix = std::array<std::array<double,gts>,alleles>;
// using GT_Matrix = std::array<std::array<int,Genotype::enumerate_gts(ploidy).size()>,Genotype::alleles.size()>;

class Popstatem{
public:
	typedef std::tuple<double, std::map<char, double>, double, double> theta_t; //theta, pi, refweight, epsilon	
protected:
	Pileupdata plp;
	theta_t theta;
	EM<double,std::map<char,double>,double,double> em;
	int ploidy;
	GT_Matrix m;
	std::vector<Genotype> possible_gts;
public:
	Popstatem(Pileupdata p, int ploidy);
	Popstatem(Pileupdata p);
	Popstatem(std::string samfile, std::string refname);
	Popstatem(std::string samfile, std::string refname, int ploidy);
	theta_t start(double stop);
	double q_function(theta_t theta);
	theta_t m_function(theta_t theta);
	void load_matrix(GT_Matrix &m, std::vector<char> x, char ref);
	void apply_over_gt(std::function<void (int, int, std::vector<Genotype>::iterator)> f);
	double dq_dtheta(double th);
	double ddq_dtheta(double th);
	double dq_dw(double w);
	double ddq_dw(double w);
	double dq_dpi(char a, double pi);
	double ddq_dpi(char a, double pi);
	static double allele_alpha(char allele, char ref, double ref_weight, double theta, std::map<char,double> pi);
	static double allele_alpha(char allele, char ref, double ref_weight, double theta, double pi);
	static double ref_alpha(double ref_weight, double theta);
	static double pg_x_given_theta(Genotype g, std::vector<char> x, theta_t theta); //not log space
	static double pdata_given_theta(std::vector<char> x, theta_t theta, std::vector<Genotype> possible_gts);

};


#endif
