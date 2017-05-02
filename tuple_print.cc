#include "tuple_print.h"

template<typename T0, typename...T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T0, T...> t){
	os << '(' << std::get<0>(t); // << quote << std::get<0>(t) << quote;
	TuplePrinter<1>::print(os,t);
	return os << ')';
}

std::ostream& operator<<(std::ostream& os, const std::tuple<>){
	return os << "()";
}
