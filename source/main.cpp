/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include "dfs.hpp"
#include "triconnectivity.hpp"
#include "utilities.hpp"
#include "testing.hpp"
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

	ugraph g = *schmidt_sample_graph();
	to_file(g,"blubb");
	schmidt_is_triconnected(g);
//
//	if (argc > 1)
//	{	std::fstream f(argv[1], std::ios::in);
//		cout << "Testing separation pairs " << std::endl;
//		test_separation_pairs(f);
//		f.close();
//	} else {
//		test_separation_pairs(std::cin);
//	}
//
//	for(int i=0; i<10; i++) {
//		auto_ptr<ugraph> g = triconnected_graph(50000);
//		node s1,s2;
//		std::cout << "go" << std::endl;
//		hopcroft_tarjan_is_triconnected_nc(*g,s1,s2);
//	}



    return 0;
}
