#include <functional>
#include <tuple>
#include <stdexcept>
#include <string>
#include <iostream>
#include "em.h"

template<typename...T>
EM<T...>::EM(std::function<long double(std::tuple<T...>)> q_function, std::function<std::tuple<T...>(std::tuple<T...>)> m_function, std::tuple<T...> theta) : q_function(q_function), m_function(m_function), likelihood(0), theta(theta){
}

template<typename...T>
long double EM<T...>::likelihood_diff(long double previous, long double current){
	if (current < previous){
		throw std::logic_error("likelihood went down. previous value: " + std::to_string(previous) + ", current value: " + std::to_string(current));
	} else {
		return current - previous;
	}
}

template<typename...T>
std::tuple<T...> EM<T...>::start(long double stop){
	long double difference;
	do{
		long double currentlike = q_function(theta);
		std::clog << "Theta: " << theta << "likelihood: " << currentlike << std::endl;
		difference = likelihood_diff(likelihood, currentlike);
		likelihood = currentlike;
		theta = m_function(theta);
	} while (difference < stop);
	return theta;
}

template<typename...T>
long double EM<T...>::get_likelihood(){
	return likelihood;
}



