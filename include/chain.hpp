/*
 * chain.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CHAIN_HPP_
#define CHAIN_HPP_

#include "LEDA/core/slist.h"
#include "LEDA/graph/ugraph.h"

using namespace leda;

enum chain_type {one, two_a, two_b, three_a, three_b, unmarked};
struct chain {
	unsigned int number;
	chain* parent;
	int segment;
	slist<chain> children;
	node s;
	node t;
	chain_type type;
	edge first_edge; // for traversals of the chain
	bool is_backedge;
	bool is_marked;
	bool in_subdivision;
	void set_parent(chain* p);
	chain();
	bool operator==(chain& c);

};

#endif /* CHAIN_HPP_ */
