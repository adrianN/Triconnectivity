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
class chain {
private:
	chain* parent;
	node s;
	node t;
public:
	unsigned int number;
	slist<chain*> children;

	chain_type type;
	edge first_edge; // for traversals of the chain
	bool is_backedge;
	bool is_marked;
	bool in_subdivision;
	void set_parent(chain* p);
	chain* get_parent(void) const;
	void set_s(node s);
	node get_s(void) const;
	void set_t(node t);
	node get_t(void) const;
	chain();
	bool operator==(chain& c) const;

};

std::ostream& operator <<(std::ostream& c, const chain_type& t);
std::ostream& operator <<(std::ostream& c, chain* const& t);
unsigned int get_number(chain* const & c);
#endif /* CHAIN_HPP_ */
