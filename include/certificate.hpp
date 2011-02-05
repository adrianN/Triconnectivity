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
#include <vector>
using namespace leda;




class certificate {
private:
	ugraph& the_graph;
public:
	certificate(ugraph& graph);
	bool add_bg_path(std::vector<edge> const & edges);
	bool add_bg_path(list<edge> const & edges);

	bool add_bg_path(const chain * a_chain);
	bool verify(void) const;
};

#endif /* CERTIFICATE_HPP_ */
