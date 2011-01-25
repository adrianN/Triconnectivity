#include "chain.hpp"

void chain::set_parent(chain* p) {
	parent = p;
	p->children.append(*this);
}
chain::chain() : parent(NULL), segment(-1), s(NULL), t(NULL),  type(unmarked), is_backedge(true), is_marked(false), in_subdivision(false) {}
bool chain::operator==(chain& c) {
	return number == c.number;
}
