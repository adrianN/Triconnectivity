#include "not_triconnected_exception.hpp"

not_triconnected_exception::not_triconnected_exception(string s) :connectivity(0) {std::cout << s << std::endl;}
not_triconnected_exception::not_triconnected_exception(string s, node art_point) : connectivity(1), articulation_point(art_point) {std::cout << s << std::endl;}
not_triconnected_exception::not_triconnected_exception(string s, std::pair<node,node> sep_pair) : connectivity(2) {
	std::cout << s<< std::endl;
	separation_pair[0] = sep_pair.first;
	separation_pair[1] = sep_pair.second;
}
