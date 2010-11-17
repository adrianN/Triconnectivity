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
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <memory>

using namespace leda;
using std::cout;
using std::endl;
using std::ostream;
using std::auto_ptr;

auto_ptr<ugraph> sample_graph(void) {
    auto_ptr<ugraph> g(new ugraph());
    node nodes[13];
    for(unsigned int i=0; i<13; i++) {
        nodes[i] = g->new_node();
    }
    int edges[13][6] = {
        {1,3,7,11,12,-1},
        {12,2,-1,-1,-1,-1},
        {12,3,-1,-1,-1,-1},
        {4,5,6,-1,-1,-1},
        {5,6,7,-1,-1,-1},
        {6,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {10,8,11,-1,-1,-1},
        {9,10,11,-1,-1,-1},
        {10,11,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<13; i++) {
        for(unsigned int j=0; j<6; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

auto_ptr<ugraph> k4(void) {
	const unsigned int num_nodes = 4;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,3,2},
    		{2,3,-1},
    		{3,-1,-1},
    		{-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}
/* K4 - one edge */
auto_ptr<ugraph> test_eins(void) {
	const unsigned int num_nodes = 4;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,3,-1},
    		{2,3,-1},
    		{3,-1,-1},
    		{-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

auto_ptr<ugraph> test_zwei(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 5;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,2,3,-1,-1},
    		{2,3,4,5,-1},
    		{3,4,5,-1,-1},
    		{-1,-1,-1,-1,-1},
    		{5,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

auto_ptr<ugraph> test_drei(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 5;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,2,3,4,-1},
    		{2,3,4,5,-1},
    		{3,4,5,-1,-1},
    		{-1,-1,-1,-1,-1},
    		{5,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

void to_file(const ugraph& g, const char* name) {
    std::ostringstream s;
    s << "./graph" << name << ".dot";
    std::fstream f(s.str().c_str(),std::ios_base::out);
    simple_to_dot(g, f);
}

void plantri_test() {
    ugraph g;
    node s1,s2;
    unsigned int graph_number = 0;
    while(!std::cin.eof()) {

    	graph_number++;
    	std::cin >> g;
    	if (!hopcroft_tarjan_is_triconnected_nc(g,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./graph" << graph_number << ".dot";
			std::fstream dot(s.str().c_str(), std::ios::out);
			simple_to_dot(g,dot);
			dot.close();
    	}
    	if (graph_number % 10000 == 0)
    		cout << '.';
    	if (graph_number % 100000 == 0)
    		cout << '\n';
    	assert(g.number_of_nodes() > 0);
    }

    cout << graph_number << endl;
}

void test_separation_pairs(std::istream& cin) {

    ugraph one,two;
    node s1=NULL,s2=NULL;
    node found1=NULL, found2=NULL;
    unsigned int graph_number = 0;
    while(!cin.eof()) {
    	s1=NULL;
    	s2=NULL;
    	found1=NULL;
    	found2=NULL;
    	graph_number++;
    	cin >> one;
    	if (cin.eof()) return;
    	cin >> two;
    	cout << "-------" << endl;
    	glue_graphs(one,two,s1,s2);


    	to_file(one,"initial");

    	const bool triconnected = hopcroft_tarjan_is_triconnected_nc(one,found1,found2);
    	bool correct_pair = false;
    	if (!triconnected)
    		correct_pair = !((s1->id()==found1->id() && s2->id()==found2->id()) || (s1->id() == found2->id() && s2->id() == found1->id()));

    	if ( triconnected || correct_pair) {
    		cout << "drama "<< graph_number << " " << triconnected << " " << correct_pair << endl;
    		if (found1!=NULL && found2!=NULL)
				cout << " found " << found1->id() << " " << found2->id()  << endl;
    		else
    			cout << " no pair found" << endl;
			std::ostringstream s;
			s << "./graph" << graph_number << ".dot";
			std::fstream dot(s.str().c_str(), std::ios::out);
			simple_to_dot(one,dot);
			dot.close();
			break;
    	}

    	if (graph_number % 10000 == 0)
    		cout << '.';
    	if (graph_number % 100000 == 0)
    		cout << '\n';
    }


}

// not triconnected
auto_ptr<ugraph> test_four(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 4;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,4,5,-1},
    		{2,3,5,-1},
    		{3,4,-1,-1},
    		{4,-1,-1,-1},
    		{5,-1,-1,-1},
    		{-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

int main(void) {
//	std::fstream f("./out.txt", std::ios::in);
//	test_separation_pairs(f);

	auto_ptr<ugraph> test = test_four();
	node s1=NULL,s2=NULL;
	bool tc = hopcroft_tarjan_is_triconnected_nc(*test,s1,s2);
	cout << "tc " << tc << endl;
	if (!tc) {
		cout << s1->id() << endl;
		cout << s2->id() << endl;
	}

    return 0;
}
