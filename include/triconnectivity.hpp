/*
 * triconnectivity.hpp
 *
 *  Created on: Oct 25, 2010
 *      Author: adrian
 */

#ifndef TRICONNECTIVITY_HPP_
#define TRICONNECTIVITY_HPP_

#include "LEDA/graph/ugraph.h"
#include "LEDA/core/tuple.h"
#include "LEDA/core/list.h"
#include <memory>


bool naive_is_triconnected(leda::ugraph&);
std::auto_ptr<leda::list<leda::two_tuple<leda::node,leda::node> > > naive_separation_pairs(leda::ugraph& g);

#endif /* TRICONNECTIVITY_HPP_ */
