#include "not_triconnected_exception.hpp"

not_triconnected_exception::not_triconnected_exception() :connectivity(0) {}
not_triconnected_exception::not_triconnected_exception(node art_point) : connectivity(1), articulation_point(art_point) {}
not_triconnected_exception::not_triconnected_exception(std::pair<node,node> sep_pair) : connectivity(2) {
	separation_pair[0] = sep_pair.first;
	separation_pair[1] = sep_pair.second;
}
