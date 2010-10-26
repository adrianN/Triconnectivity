#include "triconnectivity.hpp"
#include "LEDA/graph/node_array.h"
#include "LEDA/graph/edge_array.h"

using namespace leda;
/*
 * Finding separation pairs as described in "Dividing a graph into triconnected components"
 */

class palm_tree {
public:

	palm_tree(const ugraph& g) :
		the_graph(g),
		number(g,-1),
		flag(g,false),
		is_frond(g,false),
		lowpoint_one(new unsigned int[g.number_of_nodes()]),
		lowpoint_two(new unsigned int[g.number_of_nodes()]),
		number_descendants(new unsigned int[g.number_of_nodes()]),
		father(new unsigned int[g.number_of_nodes()])
	{
		n=0;
		dfs(the_graph.first_node(), NULL);
	};

	~palm_tree() {
		delete[] lowpoint_one;
		delete[] lowpoint_two;
		delete[] number_descendants;
		delete[] father;
	};

private:
	void dfs(const node& v,  node const * const u) {

		number[v] = n++;

		//dummy statement

		edge wv;
		forall_adj_edges(wv,v) {
			node w = target(wv);
			assert(v==source(wv));
			if (number[w] == -1 || (number[w]==0 && w!=the_graph.first_node())) { //latter part is fix for bug in LEDA
				dfs(w, &v);

				//dummy statement

			} else if (number[w] < number[v] && (u == NULL || (w!=*u) || flag[v]) ) {
				is_frond[wv] = true;

				//dummy statement
			}

			if (u!=NULL && w==*u) flag[v] = true;
		}
	};

	unsigned int n;
	const ugraph& the_graph;
	node_array<int> number; //stores the step in which a node is reached during DFS
	node_array<bool> flag; //becomes true when the entry in adj(v) corresponding to tree arc (u,v) is examined
	edge_array<bool> is_frond; //true if edge is a frond.

	unsigned int* const lowpoint_one; //stores for each node the first lowpoint
	unsigned int* const lowpoint_two;
	unsigned int* const number_descendants;
	unsigned int* const father;

};
