#ifndef __MEEP_TUPLE_PRINT_INCLUDED__
#define __MEEP_TUPLE_PRINT_INCLUDED__

#include <tuple>
#include <iostream>
#include <string>
#include <type_traits>

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
std::ostream& operator<<(std::ostream& os, const std::tuple<T0, T...> t);

std::ostream& operator<<(std::ostream& os, const std::tuple<>);


#endif
