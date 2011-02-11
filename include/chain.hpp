/*
 * chain.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CHAIN_HPP_
#define CHAIN_HPP_

#include "LEDA/core/list.h"
#include "LEDA/graph/ugraph.h"

using namespace leda;

enum chain_type {one, two_a, two_b, three_a, three_b, unmarked};
class chain {
private:
	chain* parent;
	node s;
	node t;
	list<chain*>* delete_from_list;
	list<chain*>::item delete_from_item;
	bool in_subdivision;
	chain_type the_type;
public:

	unsigned int number;
	list<chain*> children12;
	list<chain*> type3;


	edge first_edge; // for traversals of the chain
	bool is_backedge;
	bool is_marked;
	void set_parent(chain* p);
	chain* get_parent(void) const;
	void set_s(node s);
	node get_s(void) const;
	void set_t(node t);
	node get_t(void) const;
	bool is_in_subdivision(void) const;
	void add_to_subdivision(void);
	void add_to_type3(chain* c);
	void set_type(chain_type t);
	chain_type get_type(void) const;
	chain();
	bool operator==(chain& c) const;

};

std::ostream& operator <<(std::ostream& c, const chain_type& t);
std::ostream& operator <<(std::ostream& c, chain* const& t);
unsigned int get_number(chain* const & c);
#endif /* CHAIN_HPP_ */
