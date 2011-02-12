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
	fstream f;

	if (argc>1)
		{
			cout << "reading file " << argv[1] << std::endl;

			f.open(argv[1], std::ios::in);
		}
	else {
		cout << "reading file ./9.txt" << std::endl;
		f.open("./9.txt", std::ios::in);
	}

//	ugraph g;
////	g.read("./the_graph.leda");
//	f >> g;
//	g.permute_edges();
//	node s3 = NULL,s4 = NULL;
//	hopcroft_tarjan_is_triconnected_nc(g,s3,s4);
//	node n;
//	forall_nodes(n,g) {
//		node s1 = NULL, s2 = NULL;
//
//		std::cout << "Starting node " << n->id() << std::endl;
//	schmidt_is_triconnected(g,s1,s2,n);
//		if (s1!=NULL && s2 != NULL) {
//			std::cout <<s1->id() << " " << s2->id() << std::endl;
//			std::cout << s3->id() << " " << s4->id() << std::endl;
//			if (!((s1->id() == s3->id() && s2->id() == s4->id()) || (s1->id() == s4->id() && s2->id() == s3->id()))) {
//				exit(-1);
//			}
//		} else {
//			std::cout << "none found" << std::endl;
//
//			std::cout << (s3!=NULL ? s3->id() : -1) << " " << (s4!=NULL? s4->id():-1) << std::endl;
//			exit(-1);
//		}
//	}

	plantri_schmidt_test(f);
	//schmidt_test_separation_pairs(f);
    return 0;
}
