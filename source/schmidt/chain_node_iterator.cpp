#include "chain_node_iterator.hpp"

chain_node_iterator::chain_node_iterator(const chain* c, const schmidt_triconnectivity* t) : the_chain(c), decomposition(t), current_node(c->s), done(false) {
	assert(c!=NULL);
	assert(t!=NULL);
}

bool chain_node_iterator::has_next() const {
	return !done;
}

node chain_node_iterator::next() {
	assert(has_next());
	if (current_node==the_chain->t) done=true;

	node to_return = current_node;
	if (current_node == the_chain->s) {
		assert(the_chain->first_edge!=NULL);
		current_node = opposite(the_chain->s, the_chain->first_edge);
	} else {
		current_node = decomposition->parent_node(current_node);
	}

	return to_return;
}

