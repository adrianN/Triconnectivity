#include "chain.hpp"

chain::chain() : parent(NULL), s(NULL), t(NULL), delete_from_list(NULL),delete_from_item(NULL), in_subdivision(false), the_type(unmarked), first_edge(NULL), is_backedge(true), is_marked(false) {}
bool chain::operator==(chain& c) const {
	return number == c.number;
}

unsigned int get_number(chain* const & c) { assert(c!=NULL); return c->number; }

std::ostream& operator <<(std::ostream& c, chain* const& t) {
	if (t==NULL)
		c << "null";
	else
		c << t->number;
	return c;
}


std::ostream& operator <<(std::ostream& c, const chain_type& t) {
	switch(t) {
	case one: c<<"1"; break;
	case two_a: c<<"2a"; break;
	case two_b: c<<"2b"; break;
	case three_a: c<<"3a"; break;
	case three_b: c<<"3b"; break;
	case unmarked: c<<"X"; break;
	};
	return c;
}
