/*
 * not_triconnected_exception.hpp
 *
 *  Created on: Feb 8, 2011
 *      Author: adrian
 */

#ifndef NOT_TRICONNECTED_EXCEPTION_HPP_
#define NOT_TRICONNECTED_EXCEPTION_HPP_

#include "LEDA/graph/ugraph.h"
#include <string>

using namespace leda;

class not_triconnected_exception {
public:
	not_triconnected_exception(std::string s);
	not_triconnected_exception(std::string s,node art_point);
	not_triconnected_exception(std::string s,node sep1, node sep2);
	unsigned int connectivity;
	union {
		node articulation_point;
		node separation_pair[2];
	};
	std::string message;
};

#endif /* NOT_TRICONNECTED_EXCEPTION_HPP_ */
