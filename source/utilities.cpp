#include "triconnectivity.hpp"
#include <ostream>
#include <vector>

using namespace leda;
using std::endl;
using std::istream;

void to_dot(const ugraph& g, std::ostream& out) {

    out << "digraph G {" << endl;
    node n;
    forall_nodes(n,g) {

        out << "\tnode" << n->id() << " [label = \"" << n->id() << "\"]" << endl;

        node u;
        unsigned int i=0;
        forall_adj_nodes(u,n) {
            out << "\tnode" << n->id() << " -> " << "node" << u->id() << " [label = \"" << ++i << "\"]" << endl;
        }
    }

    out << "}" << endl;
}

void simple_to_dot(const ugraph& g, std::ostream& out) {
    out << "graph G {" << endl;
    {   edge e;
        forall_edges(e,g) {
            out << "\tnode" << source(e)->id() << " -- " << "node" << target(e)->id() << endl;
        }
    }
    {   node n;
        forall_nodes(n,g) {
            out << "\tnode" << n->id() << " [label=\""<<n->id()<<"\"]" << endl;
        }
    }
    out << "}" << endl;
}

istream& operator>>(istream& input, ugraph& graph) {
//	std::cout << "Graph input " << input.tellg() << std::endl;
	if (input.tellg() == std::ios::beg) {
		input.ignore(15); // the file starts with >>planar_code<<
	}
	graph.clear();
	unsigned int num_nodes;
	assert(input.good());
	assert(!input.eof());
	num_nodes = input.get();
	if (num_nodes == '>' && input.peek() == '>') {
		input.ignore(14);
		num_nodes = input.get();
		assert(num_nodes != '<');
	}

//	std::cout << "\t number of nodes: " << num_nodes << std::endl;
	std::vector<node> nodes(num_nodes);

	for(unsigned int i=0; i<nodes.size(); i++) {
		nodes[i] = graph.new_node();
	}


	for(unsigned int i=0; i<num_nodes; i++) {
		unsigned int adj_node;
		adj_node = input.get();
		while (adj_node!=0) {
//			std::cout<< "\tEdge from " << i << " to " << (int)adj_node-1 << std::endl;
			if (adj_node-1>i)
				graph.new_edge(nodes[i],nodes[adj_node-1]);
			adj_node = input.get();
			assert(input.good());
		}
	}
	input.peek();
	return input;
}

node merge_nodes(ugraph& g, node one, node two) {
	node new_node = g.new_node();
	node_array<bool> edge_to(g,false);
	edge e;
	forall_adj_edges(e,one) {
		node opp = opposite(e,one);
		if (!edge_to[opp]) {
			g.new_edge(new_node,opp);
			edge_to[opp] = true;
		}
	}
	forall_adj_edges(e,two) {
		node opp = opposite(e,two);
		if (!edge_to[opp]) {
			g.new_edge(new_node,opp);
			edge_to[opp] = true;
		}
	}
	g.del_node(one);
	g.del_node(two);

	return new_node;
}

/* Takes two graphs and indentifies two nodes from each to form a new graph with exactly one separation pair. The result is stored in one */
void glue_graphs(ugraph& one, ugraph& two) {
	const unsigned int number_two_nodes = two.number_of_nodes();

	const node point_one = one.choose_node();
	node point_two;
	do {
		point_two = one.choose_node();
	} while(point_one == point_two);

	node* const new_nodes = new node[number_two_nodes];
	for(unsigned int i = 0; i< number_two_nodes; i++) {
		new_nodes[i] = one.new_node();
	}

	edge e;
	forall_edges(e,two) {
		const unsigned int u = source(e)->id();
		const unsigned int v = target(e)->id();
		one.new_edge(new_nodes[u],new_nodes[v]);
	}

	random_source src_two(0,number_two_nodes-1);

	unsigned int other_graph_point_one, other_graph_point_two;
	src_two >> other_graph_point_one;
	do {
		src_two >> other_graph_point_two;
	} while(other_graph_point_one == other_graph_point_two);

	node merge_one = merge_nodes(one,point_one, new_nodes[other_graph_point_one]);
	node merge_two = merge_nodes(one,point_two, new_nodes[other_graph_point_two]);

	std::cout << "Merge nodes: " << merge_one->id() << " " << merge_two->id() << " from " << std::endl;
	std::cout << point_one->id() << " " << new_nodes[other_graph_point_one]->id() << std::endl;
	std::cout << point_two->id() << " " << new_nodes[other_graph_point_two]->id() << std::endl;

	delete[] new_nodes;
}



