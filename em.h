#ifndef __MEEP_EM_INCLUDED__
#define __MEEP_EM_INCLUDED__

#include <tuple>
#include <vector>

//see https://www2.ee.washington.edu/techsite/papers/documents/UWEETR-2010-0002.pdf for a tutorial on EM
template<typename... T>
class EM{
protected:
	long double likelihood;
	std::tuple<T...> theta;
	std::function<long double(std::tuple<T...> theta)> q_function; //returns expected value of log likelihood function
	std::function<std::tuple<T...>(std::vector<T...> theta)> m_function; //returns theta that maximizes Q
	EM(std::function<long double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta); //initialize with guess for theta
	long double likelihood_diff(long double, long double);
public:
	std::tuple<T...> start(long double stop); //start the EM. return theta.
	long double get_likelihood();
};

#endif