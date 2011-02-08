#include "chain_edge_iterator.hpp"

chain_edge_iterator::chain_edge_iterator(const chain* c, const schmidt_triconnectivity* t) : the_chain(c), decomposition(t), current_edge(c->first_edge) {
	assert(c!=NULL);
	assert(t!=NULL);
}

edge chain_edge_iterator::next(void) {
	if (current_edge == NULL) return NULL;
	edge to_return = current_edge;
	const node endnode = the_chain->get_t();

	if (target(current_edge) == endnode) {
		current_edge = NULL;
		return to_return;
	}

	current_edge = decomposition->parent[target(current_edge)];

	return to_return;
}


