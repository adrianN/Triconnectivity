/*
 * utilities.hpp
 *
 *  Created on: Oct 25, 2010
 *      Author: adrian
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include "LEDA/graph/ugraph.h"
#include <ostream>
#include <istream>
#include <memory>

void to_dot(const leda::ugraph& g, std::ostream& output);
void simple_to_dot(const leda::ugraph& g, std::ostream& output);

std::istream& operator>>(std::istream& input, leda::ugraph& graph);
leda::node merge_nodes(leda::ugraph& g, leda::node one, leda::node two);
void glue_graphs(leda::ugraph& one, leda::ugraph& two, leda::node& merge_one, leda::node& merge_two);
void write_planar_code(leda::ugraph& g, std::ostream& out);
std::auto_ptr< leda::ugraph > triconnected_graph(const unsigned int n);

#endif /* UTILITIES_HPP_ */
