#ifndef __MEEP_EM_INCLUDED__
#define __MEEP_EM_INCLUDED__

//see https://www2.ee.washington.edu/techsite/papers/documents/UWEETR-2010-0002.pdf for a tutorial on EM
class EM{
protected:
	long double likelihood;
	std::vector<long double> theta;
	std::function<long double(std::vector<long double> theta)> q_function; //returns expected value of log likelihood function
	std::function<std::vector(std::vector<long double> theta)> m_function; //returns theta that maximizes Q
	EM(std::function<long double(std::vector)> q_function, std::function<std::vector(std::vector(std::vector)> m_function, std::vector theta); //initialize with guess for theta
	long double likelihood_diff(long double, long double);
public:
	std::vector start(long double stop); //start the EM. return theta.
	long double get_likelihood();
};

#endif