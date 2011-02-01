#include "chain_node_iterator.hpp"

chain_node_iterator::chain_node_iterator(const chain* c, const schmidt_triconnectivity* t) : the_chain(c), decomposition(t), current_node(c->get_s()) {
	assert(c!=NULL);
	assert(t!=NULL);
}



node chain_node_iterator::next() {


	node to_return = current_node;

	if (current_node == NULL || current_node==the_chain->get_t()) {
		current_node = NULL;
		return to_return;
	}

	if (current_node == the_chain->get_s() && the_chain->number!=0) {
		assert(the_chain->first_edge!=NULL);
		current_node = target(the_chain->first_edge);
	} else {
		current_node = decomposition->parent_node(current_node);
	}

	return to_return;
}

