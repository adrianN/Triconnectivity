/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include <iostream>
#include <ostream>
#include <fstream>
#include <cassert>

using namespace leda;
using std::cout;
using std::endl;
using std::ostream;

/* Does a DFS to determine whether the graph is connected or not */
bool is_connected(const ugraph& g) {
	assert(g.number_of_nodes()>=0);
	unsigned int number_of_nodes = g.number_of_nodes();
	node* stack = new node[number_of_nodes];
	node* base = stack;
	*stack = g.first_node();

	node_array<bool> visited(g,false);
	visited[*stack] = true;

	unsigned int number_seen=1;

	while(stack>=base) {
		/*const*/ node current = *stack;
		stack--;
		node next_node;
		forall_adj_nodes(next_node,current) {
			if (!visited[next_node]) {
				stack++;
				*stack = next_node;
				visited[next_node] = true;
				number_seen++;
			}
		}
	}

	delete base;
	return number_seen==number_of_nodes;
}

bool naive_is_triconnected(const ugraph& g) {
	/* A graph is triconnected if there is no pair of nodes u,v s.t.
	 * the graph is no longer connected if u and v are removed. i.e.
	 * there exists a pair of nodes x,y s.t. all paths between x,y pass through u or v.
	 */
	return true;
}

void to_dot(const ugraph& g, ostream& out) {
	node_array<int> seen(g,0); //initialisation doesn't work
	edge e;
	unsigned int number=1;

	out << "graph G {" << endl;

	forall_edges(e,g) {
		node src = source(e);
		node trg = target(e);
		if (seen[src] == 0) {
			seen[src] = number++;
		}
		if (seen[trg] == 0) {
			seen[trg] = number++;
		}

		out << "\tnode" << seen[src] << " -- " << "node" << seen[trg] << endl;
	}

	node n;
	forall_nodes(n,g) {
		out << "\tnode"<<seen[n]<<" [label= "<<seen[n]<<"]" << endl;
	}

	out << "}" << endl;
}

int main() {
	ugraph g;
	const double p = 0.5;
	const unsigned int n = 1000;
	const unsigned int m = (int)(((double)(n*(n-1))*p)/2.0);
	float total = used_time();
	for(int i=0; i<10; i++) {
		random_graph(g,n,m,false,true,false);
		//std::fstream f("./leDot.dot",std::ios_base::out);
		//to_dot(g,f);
//		float ut = used_time();
		cout << is_connected(g) << endl;
//		cout << used_time()-ut;
	}
	cout << "total " << used_time()-total << endl;
	return 0;
}
