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
#include <sstream>
#include <string>
#include <cassert>

using namespace leda;
using std::cout;
using std::endl;
using std::ostream;

void to_dot(const ugraph&,ostream&);

/* Does a DFS to determine whether the graph is connected or not */
bool is_connected(const ugraph& g) {
	assert(g.number_of_nodes()>=0);
	unsigned int number_of_nodes = g.number_of_nodes();
	node* stack = new node[number_of_nodes];
	node* base = stack;
	*stack = g.first_node();
	assert(!g.is_hidden(*stack));

	node_array<bool> visited(g,false);
	assert(visited[*stack]==false);
	visited[*stack] = true;

	unsigned int number_seen=1;

	while(stack>=base) {
		/*const*/ node current = *stack;
		stack--;
		g.print_node(current);
		cout << " ("<< g.degree(current)<<") ";
		node next_node;
		forall_adj_nodes(next_node,current) {
			assert(!g.is_hidden(next_node));
			cout << ',';
			if (!visited[next_node]) {
				cout << '.';
				stack++;
				*stack = next_node;
				visited[next_node] = true;
				number_seen++;
			}
		}
	}

	cout << endl;

	delete base;
	cout << "reachable " << number_seen << endl;
	return number_seen==number_of_nodes;
}

bool naive_is_triconnected(ugraph& g) {
	/* A graph is triconnected if there is no pair of nodes u,v s.t.
	 * the graph is no longer connected if u and v are removed. i.e.
	 * there exists a pair of nodes x,y s.t. all paths between x,y pass through u or v.
	 */
	node u;
	forall_nodes(u,g) {
		int degree = g.degree(u);
		list<edge> adj_edges = g.adj_edges(u);
		g.hide_node(u);

		node v;
		forall_nodes(v,g) {
			int degree = g.degree(v);
			list<edge> adj_edges = g.adj_edges(v);
			g.hide_node(v);

			if (!is_connected(g)) { //u,v is a separation pair
				g.print_node(u);
				g.print_node(v);
				cout << endl;
//				std::fstream f("./blubb.dot",std::ios_base::out);
//				to_dot(g,f);
//				f.close();
				g.restore_all_nodes();
				g.restore_all_edges();
				return false;
			}
			g.restore_node(v);
			g.restore_edges(adj_edges);
			assert(!g.is_hidden(v));
			assert(degree==g.degree(v));

		}
		g.restore_node(u);
		g.restore_edges(adj_edges);
		assert(degree==g.degree(u));

	}
	return true;
}

void to_dot(const ugraph& g, ostream& output) {
	edge e;
	std::ostringstream out;
	out << "graph G {" << endl;

	forall_edges(e,g) {
		node src = source(e);
		node trg = target(e);

		out << "\tnode";
		g.print_node(src,out);
		out << " -- node";
		g.print_node(trg,out);
		out << endl;
	}

	out << "}" << endl;
	std::string s = out.str();
	for(unsigned int i=0; i<s.length(); i++) {
		if (s[i]=='[' || s[i] == ']') {
			s[i] = '_';
		}
	}

	output << s;
}

int main(void) {
	ugraph g;
	const double p = 0.5;
	const unsigned int n = 10;
	const unsigned int m = (int)(((double)(n*(n-1))*p)/2.0);
	cout << "n, m " << n << ", "<<m<<endl;
	float total = used_time();
	for(int i=0; i<1; i++) {
		random_graph(g,n,m,false,true,false);

		std::ostringstream s;
		s << "./graph" << i << ".dot";
		std::fstream f(s.str().c_str(),std::ios_base::out);
		to_dot(g, f);
		cout << "==" <<i << "==" <<endl;
		cout << is_connected(g) << endl;
		cout << naive_is_triconnected(g) << endl;


	}
	cout << "total " << used_time()-total << endl;
	return 0;
}
