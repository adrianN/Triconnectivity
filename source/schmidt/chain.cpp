#include "chain.hpp"

void chain::set_s(node s) {
	assert(this->s==NULL);
	this->s = s;
}
node chain::get_s(void) const {
	assert(s!=NULL);
	return s;
}

void chain::set_t(node t) {
	assert(this->t==NULL);
	this->t=t;
}

node chain::get_t(void) const {
	assert(t!=NULL);
	return t;
}

chain* chain::get_parent(void) const {
	return parent;
}

void chain::set_parent(chain* p) {
	assert(parent==NULL);
	parent = p;
	if (p!=NULL)
		p->children.append(this);
	else
		assert(number==0);
}
chain::chain() : parent(NULL), s(NULL), t(NULL), type(unmarked), first_edge(NULL), is_backedge(true), is_marked(false), in_subdivision(false) {}
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
