/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include <LEDA/core/slist.h>
#include "dfs.hpp"
#include "triconnectivity.hpp"
#include "utilities.hpp"
#include "testing.hpp"
#include "InterlockingIntervals.hpp"
#include <iostream>

using namespace leda;
using namespace std;

void to_file(const ugraph& g, const char* name) {
    std::ostringstream s;
    s << "./graph" << name << ".dot";
    std::fstream f(s.str().c_str(),std::ios_base::out);
    simple_to_dot(g, f);
}






int main(int argc, char* argv[]) {
	plantri_schmidt_test(cin);
    return 0;
}
