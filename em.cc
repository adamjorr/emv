EM::EM(std::function<long double(std::vector)> q_function, std::function<std::vector(std::vector(std::vector)> m_function, std::vector theta) : q_function(q_function), m_function(m_function), likelihood(0), theta(theta){
}

long double EM::likelihood_diff(long double previous, long double current){
	if (current < previous){
		throw std::runtime_error("likelihood went down. previous value: " + std::to_string(previous) ", current value: " + std::to_string(current));
	} else {
		return current - previous;
	}
}

std::vector EM::start(long double stop){
	do{
		long double currentlike = q_function(theta);
		long double difference = likelihood_diff(likelihood, currentlike);
		likelihood = currentlike;
		theta = m_function(theta);
	} while (difference < stop);
	return theta;
}

long double EM::get_likelihood(){
	return likelihood;
}


