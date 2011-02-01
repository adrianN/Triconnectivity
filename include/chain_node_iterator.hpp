/*
 * chain_node_iterator.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CHAIN_NODE_ITERATOR_HPP_
#define CHAIN_NODE_ITERATOR_HPP_

#include "chain.hpp"
#include "schmidt.hpp"

class chain_node_iterator {
private:
	const chain* the_chain;
	const schmidt_triconnectivity* decomposition;
	node current_node;
public:
	chain_node_iterator(const chain* c, const schmidt_triconnectivity* t);
	node next();
};


#endif /* CHAIN_NODE_ITERATOR_HPP_ */
