#include "tuple_print.h"

std::ostream& operator<<(std::ostream& os, const std::tuple<>){
	return os << "()";
}
