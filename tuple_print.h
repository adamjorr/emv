#include <tuple>
#include <iostream>
#include <string>

//check out http://en.cppreference.com/w/cpp/utility/tuple/tuple_cat for another example
// also http://stackoverflow.com/questions/6245735/pretty-print-stdtuple
template<typename Tuple, std::size_t i>
struct TuplePrinter{
	operator<<(std::ostream& os, const Tuple& t){

	}
};







std::ostream& operator<<(std::ostream& os, const std::tuple<>&){
	return os << "()";
}

template<typename T0, typename ...T>
std::ostream& print_tuple(std::ostream& os, const std::tuple<T0, T...>& t){
	os << std::get<0>(t) << " , " << print_tuple(os,)
}

template<typename T0, typename ...T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T0, T...>& t){
	os << '(' << std::get<0>(t) ;

}




