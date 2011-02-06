/*
 * chain_edge_iterator.hpp
 *
 *  Created on: Feb 6, 2011
 *      Author: adrian
 */

#ifndef CHAIN_EDGE_ITERATOR_HPP_
#define CHAIN_EDGE_ITERATOR_HPP_

#include "chain.hpp"
#include "schmidt.hpp"



class chain_edge_iterator {
private:
	const chain* the_chain;
	const schmidt_triconnectivity* decomposition;
	edge current_edge;
public:
	chain_edge_iterator(const chain* c, const schmidt_triconnectivity* t);
	edge next();
};


#endif /* CHAIN_EDGE_ITERATOR_HPP_ */
