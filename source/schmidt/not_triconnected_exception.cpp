#include "not_triconnected_exception.hpp"

not_triconnected_exception::not_triconnected_exception(string s) :connectivity(0), message(s) {}
not_triconnected_exception::not_triconnected_exception(string s, node art_point) : connectivity(1), articulation_point(art_point), message(s) {

}
not_triconnected_exception::not_triconnected_exception(string s, node sep1, node sep2) : connectivity(2), message(s) {
	separation_pair[0] = sep1;
	separation_pair[1] = sep2;
}
