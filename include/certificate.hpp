/*
 * certificate.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CERTIFICATE_HPP_
#define CERTIFICATE_HPP_

#include "LEDA/graph/ugraph.h"
#include "chain.hpp"
#include "construction_sequence.hpp"

using namespace leda;

union cert_data {
	unsigned int number_of_vertices;
	two_tuple<node,node>* sep_pair;
	chain* decomposition;
	construction_sequence sequence;
};

struct certificate {
private:
	ugraph& the_graph;
	unsigned int connectivity;
	cert_data* content;
public:
	certificate(ugraph& graph, unsigned int k, cert_data* c);
	bool verify();
};

#endif /* CERTIFICATE_HPP_ */
