/*
 * dfs.hpp
 *
 *  Created on: Oct 25, 2010
 *      Author: adrian
 */

#ifndef DFS_HPP_
#define DFS_HPP_

#include "LEDA/graph/ugraph.h"

bool is_connected(const leda::ugraph&);
bool dfs_order(const leda::ugraph& g, const leda::node startnode, leda::node ordered[]);


#endif /* DFS_HPP_ */
