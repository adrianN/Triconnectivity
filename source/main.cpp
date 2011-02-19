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


bool check_graph(ugraph u) {

	//simple_to_dot(u,f);
	//std::cout << "n " << u.number_of_nodes() << " m " << u.number_of_edges() << std::endl;
	node s1=NULL, s2=NULL;
	node s3=NULL, s4=NULL;
	bool s =  schmidt_is_triconnected(u,s1,s2);
	bool ht = hopcroft_tarjan_is_triconnected(u,s3,s4);
	bool all = s == ht;
	if (!all) {
		bool naive = naive_is_triconnected(u);

		to_file(u,"the_graph");
		fstream f("./the_graph.tri", ios::out);
		write_planar_code(u,f);

		std::cout << "ht " << ht << " s " << s << " naive " << naive << std::endl;
		std::cout << "Sep ";
		if (s1!=NULL) {
			u.hide_node(s1);
			std::cout << s1->id() << " ";
		} else {
			std::cout << "null ";
		}
		if (s2!=NULL) {
			u.hide_node(s2);
			std::cout << s2->id() << " ";
		} else {
			std::cout << "null ";
		}
		std::cout << std::endl;
		to_file(u,"hidden");
		assert(false);
	}

//	std::cout << "ht " << ht << " s " << s << std::endl;
	return s;
}

bool check_random_graph(int n, int m) {
	graph g;
	random_simple_undirected_graph(g,n,m);
	ugraph u(g);
	return check_graph(u);
}

void f() {
	//	fstream f;
	//
	//	if (argc>1)
	//		{
	//			cout << "reading file " << argv[1] << std::endl;
	//
	//			f.open(argv[1], std::ios::in);
	//		}
	//	else {
	//		cout << "reading file ./9.txt" << std::endl;
	//		f.open("./9.txt", std::ios::in);
	//	}
	//
	////	ugraph g;
	//////	g.read("./the_graph.leda");
	////	f >> g;
	////	g.permute_edges();
	////	node s3 = NULL,s4 = NULL;
	////	hopcroft_tarjan_is_triconnected_nc(g,s3,s4);
	////	node n;
	////	forall_nodes(n,g) {
	////		node s1 = NULL, s2 = NULL;
	////
	////		std::cout << "Starting node " << n->id() << std::endl;
	////	schmidt_is_triconnected(g,s1,s2,n);
	////		if (s1!=NULL && s2 != NULL) {
	////			std::cout <<s1->id() << " " << s2->id() << std::endl;
	////			std::cout << s3->id() << " " << s4->id() << std::endl;
	////			if (!((s1->id() == s3->id() && s2->id() == s4->id()) || (s1->id() == s4->id() && s2->id() == s3->id()))) {
	////				exit(-1);
	////			}
	////		} else {
	////			std::cout << "none found" << std::endl;
	////
	////			std::cout << (s3!=NULL ? s3->id() : -1) << " " << (s4!=NULL? s4->id():-1) << std::endl;
	////			exit(-1);
	////		}
	////	}
	//
	////	plantri_test(f);
	//	plantri_schmidt_test(f);
}

void check_planar_code(const char* filename) {
	fstream f(filename, ios::in);
	ugraph u;
	f >> u;
	node s1,s2;
	bool s =  schmidt_is_triconnected(u,s1,s2);
	bool ht = hopcroft_tarjan_is_triconnected(u,s1,s2);
	if (ht!=s) {
		std::cout << "ht " << ht << " s " << s << " naive " << naive_is_triconnected(u) << std::endl;
	} else {
		std::cout << "a-okay" << std::endl;
	}
}


void graph_statistics(int n) {
	const int num = 80;
	const int edges = (n*(n-1))/2;
	for (float f = 5; f < 25; f+=0.3) {
		int count=0;
		const float p = f/((float)n);
		const int m = (int)((float)edges*p);
		std::cout << " " << n << " " << f << " ";
		for (int i = 0; i<num; i++) {

			if (check_random_graph(n,m)) {
				count++;
			}
		}
		std::cout << (float)count/num << std::endl;
	}
}

int main(int argc, char* argv[]) {
	auto_ptr<ugraph> u = test_fc_two_children();
	to_file(*u,"simple");
	check_graph(*u);
	return 0;
}
