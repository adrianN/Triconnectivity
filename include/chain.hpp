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
	inline void set_parent(chain* p) {
	//	std::cout << "Parent of " << number << " is now " << ((p==NULL)? -1 : p->number) << std::endl;
		assert(parent==NULL);
		assert(number != 0 || p == NULL);
		parent = p;

	};
	inline chain* get_parent(void) const  {
		return parent;
	};

	inline void set_s(node s) {
		assert(this->s==NULL);
		this->s = s;
	};

	inline node get_s(void) const {
		assert(s!=NULL);
		return s;
	};

	inline void set_t(node t) {
		assert(this->t==NULL);
		this->t=t;
	};

	inline node get_t(void) const  {
		assert(t!=NULL);
		return t;
	};

	inline bool is_in_subdivision(void) const  {
		return in_subdivision;
	};

	inline void add_to_subdivision(void)  {
		in_subdivision = true;
		if (delete_from_list != NULL) {
			delete_from_list->del_item(delete_from_item);
		} else {
			assert(get_type()==unmarked);
		}
	};

	inline void add_to_type3_of(chain* c)  {
		assert(get_type() == three_a || get_type() == three_b);
		assert(delete_from_list == NULL);
		delete_from_list = &(c->type3);
		delete_from_item = delete_from_list->append(this);
	};

	inline void set_type(chain_type t)  {
		the_type = t;
		if (parent!=NULL) {
			switch(get_type()) {
				case three_a: case three_b: break;
				case unmarked: assert(false);
				default: {
					assert(delete_from_list==NULL);
					delete_from_list = &(parent->children12);
					delete_from_item = delete_from_list->append(this);
				}
			}
		}
	};

	inline chain_type get_type(void) const  {
		return the_type;
	};

	chain();
	bool operator==(chain& c) const;

};

std::ostream& operator <<(std::ostream& c, const chain_type& t);
std::ostream& operator <<(std::ostream& c, chain* const& t);
unsigned int get_number(chain* const & c);
#endif /* CHAIN_HPP_ */
