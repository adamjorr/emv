#ifndef __MEEP_TUPLE_PRINT_INCLUDED__
#define __MEEP_TUPLE_PRINT_INCLUDED__

#include <tuple>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include <map>

//check out http://en.cppreference.com/w/cpp/utility/tuple/tuple_cat for another example
// also http://stackoverflow.com/questions/6245735/pretty-print-stdtuple
template<size_t N>
struct TuplePrinter{
	template<typename...T>
	static typename std::enable_if<(N<sizeof...(T))>::type
	print(std::ostream& os, const std::tuple<T...> t){
		os << ", " << std::get<N>(t);
		TuplePrinter<N+1>::print(os,t);
	}

	template<typename...T>
	static typename std::enable_if<!(N<sizeof...(T))>::type
	print(std::ostream& os, const std::tuple<T...> t){
	}
};

template<typename T0, typename...T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T0, T...> t){
	os << '(' << std::get<0>(t); // << quote << std::get<0>(t) << quote;
	TuplePrinter<1>::print(os,t);
	return os << ')';
}

std::ostream& operator<<(std::ostream& os, const std::tuple<>);


//Vector printing
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> vec){
	if(vec.size() == 0){
		return os << "[ ]";
	}
	else{
		os << "[ " << vec[0];
		for (int i = 1; i < vec.size(); ++i){
			os << ", " << vec[i];
		}
		return os << " ]";
	}
}

//Map printing
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::map<T,U> m){
	if(m.size() == 0){
		return os << "{ }";
	}
	else{
		auto i = m.begin();
		os << "{ " << i->first << ":" << i->second;
		for(++i; i != m.end(); ++i){
			os << ", " << i->first <<  ":" << i->second;
		}
		return os << " }";
	}
}


#endif
