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

void to_dot(const leda::ugraph& g, std::ostream& output);

#endif /* UTILITIES_HPP_ */
