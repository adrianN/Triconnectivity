#include "dfs.hpp"
#include <cassert>

using namespace leda;

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
		//pop from stack
		/*const*/ node current = *stack;
		stack--;
		node next_node;

		//push unvisited neighbours
		forall_adj_nodes(next_node,current) {
			assert(!g.is_hidden(next_node));
			if (!visited[next_node]) {
				stack++;
				*stack = next_node;
				visited[next_node] = true;
				number_seen++;
			}
		}
	}

	delete[] base;
	return number_seen==number_of_nodes;
}
