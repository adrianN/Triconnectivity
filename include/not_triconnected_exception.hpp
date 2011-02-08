/*
 * not_triconnected_exception.hpp
 *
 *  Created on: Feb 8, 2011
 *      Author: adrian
 */

#ifndef NOT_TRICONNECTED_EXCEPTION_HPP_
#define NOT_TRICONNECTED_EXCEPTION_HPP_

#include "LEDA/graph/ugraph.h"

using namespace leda;

class not_triconnected_exception {
public:
	not_triconnected_exception();
	not_triconnected_exception(node art_point);
	not_triconnected_exception(std::pair<node,node> sep_pair);
	unsigned int connectivity;
	union {
		node articulation_point;
		node separation_pair[2];
	};
};

#endif /* NOT_TRICONNECTED_EXCEPTION_HPP_ */
