#ifndef __MEEP_EM_INCLUDED__
#define __MEEP_EM_INCLUDED__

#include <tuple>
#include <vector>
#include <functional>
#include <stdexcept>
#include <string>
#include <iostream>
#include "tuple_print.h"

//see https://www2.ee.washington.edu/techsite/papers/documents/UWEETR-2010-0002.pdf for a tutorial on EM
template<typename... T>
class EM{
protected:
	double likelihood;
	std::function<double(std::tuple<T...> theta)> q_function; //returns expected value of log likelihood function
	std::function<std::tuple<T...>(std::tuple<T...> theta)> m_function; //returns theta that maximizes Q	
	std::tuple<T...> theta;
	double likelihood_diff(double, double);
public:
	EM(std::function<double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta); //initialize with guess for theta
	std::tuple<T...> start(double stop); //start the EM. return theta.
	double get_likelihood();
};

//definition of template class must be in h file

template<typename...T>
EM<T...>::EM(std::function<double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta) : likelihood(0), q_function(q_function), m_function(m_function), theta(theta){
}

template<typename...T>
double EM<T...>::likelihood_diff(double previous, double current){
	if (previous == 0){
		return (current > 0 ? current : -current);
	}
	else if (current < previous){
		// throw std::logic_error("likelihood went down. previous value: " + std::to_string(previous) + ", current value: " + std::to_string(current));
		std::clog << "likelihood decreased by " << -(current - previous) << std::endl;
		return current - previous;
	} else {
		return current - previous;
	}
}

template<typename...T>
std::tuple<T...> EM<T...>::start(double stop){
	int counter = 0;
	double difference;
	do{
		double current_likelihood = q_function(theta);
		std::clog << "Theta = " << theta << "\nlikelihood = " << current_likelihood << std::endl;
		difference = likelihood_diff(likelihood, current_likelihood);
		likelihood = current_likelihood;
		theta = m_function(theta);
	// } while (difference > stop);
		counter++;
	} while (counter < 20);
	return theta;
}

template<typename...T>
double EM<T...>::get_likelihood(){
	return likelihood;
}



#endif