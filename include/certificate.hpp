/*
 * certificate.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CERTIFICATE_HPP_
#define CERTIFICATE_HPP_

#include "LEDA/graph/ugraph.h"
#include "LEDA/core/list.h"
#include "chain.hpp"
#include "schmidt.hpp"
#include <vector>
#include <utility>
using namespace leda;


class schmidt_triconnectivity; //NO IDEA why this is necessary

class certificate {
private:
	bool still_valid;
	ugraph & the_graph;
	ugraph new_graph;
	schmidt_triconnectivity* decomposition;
	node_array<int> created_by_chain;
	std::vector<list<edge>* > chains;
	edge_array<std::pair<list<edge>*, list<edge>::item>* > le_edges;
	node_array<node> orig_2_new;
	node_array<node> new_2_orig;
	std::vector<std::pair<node,node>* > endvertices;
	edge_array<bool> edge_accounted_for;

public:
	certificate(ugraph  & graph,  schmidt_triconnectivity* d);
	~certificate();
	bool add_bg_path(list<edge> const & edges) throw();

	bool add_bg_path(const chain * a_chain) throw();
	bool verify(void) throw();
};

#endif /* CERTIFICATE_HPP_ */
