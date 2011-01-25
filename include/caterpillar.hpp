/*
 * caterpillar.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CATERPILLAR_HPP_
#define CATERPILLAR_HPP_

#include "LEDA/core/slist.h"
#include "chain.hpp"

using namespace leda;

class caterpillar : public slist<chain* const> {

public:
	caterpillar(void);
	item append(chain* const & element);
	chain* parent;
};


#endif /* CATERPILLAR_HPP_ */
