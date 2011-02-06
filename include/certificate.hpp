/*
 * certificate.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef CERTIFICATE_HPP_
#define CERTIFICATE_HPP_

#include "LEDA/graph/ugraph.h"
#include "LEDA/core/list.h"
#include "chain.hpp"
#include "schmidt.hpp"
using namespace leda;


class schmidt_triconnectivity; //NO IDEA why this is necessary

class certificate {
private:
	ugraph const & the_graph;
	const schmidt_triconnectivity* decomposition;
public:
	certificate(ugraph const & graph, const schmidt_triconnectivity* d);
	bool add_bg_path(list<edge> const & edges);

	bool add_bg_path(const chain * a_chain);
	bool verify(void) const;
};

#endif /* CERTIFICATE_HPP_ */
