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
//	std::cout << "Parent of " << number << " is now " << ((p==NULL)? -1 : p->number) << std::endl;
	assert(parent==NULL);
	assert(number != 0 || p == NULL);
	parent = p;

}

bool chain::is_in_subdivision() const {
	return in_subdivision;
}

void chain::add_to_subdivision() {
	in_subdivision = true;
	if (delete_from_list != NULL) {
		delete_from_list->del_item(delete_from_item);
	} else {
		assert(get_type()==unmarked);
	}
}

void chain::add_to_type3(chain* c) {
	assert(get_type() == three_a || get_type() == three_b);
	assert(delete_from_list == NULL);
	delete_from_list = &(c->type3);
	delete_from_item = delete_from_list->append(this);
}

chain_type chain::get_type(void) const {
	return the_type;
}

void chain::set_type(chain_type t) {
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
}

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
