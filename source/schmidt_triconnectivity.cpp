#include "LEDA/graph/ugraph.h"
#include "LEDA/core/stack.h"

using namespace leda;
enum chain_type { one, two_a, two_b, three_a, three_b};

struct chain {
	chain* parent;
	node low;
	chain_type type;
};





class schmidt_triconnectivity {


	edge_array<chain> belongs_to_chain;
	edge_array<bool> is_frond;

	node_array<unsigned int> dfi;
	node_array<node> parent;
	chain* chains;
	ugraph& the_graph;

public:
	schmidt_triconnectivity(ugraph& graph) :
		belongs_to_chain(graph),
		is_frond(graph, false),
		dfi(graph,0),
		parent(graph),
		chains(new chain[graph.number_of_edges() - graph.number_of_nodes() +2]),
		the_graph(graph)
	{
		initial_dfs();
	}

	~schmidt_triconnectivity() {
		delete[] chains;
	}

	void initial_dfs() {


		    node current_node = source(the_graph.first_edge());
		    const node root = current_node;
			assert(the_graph.degree(root)>=3);

		    unsigned int number_seen=1;
		    dfi[current_node] = number_seen++;
		    parent[current_node] = NULL;
		    edge_array<bool> seen_edge(the_graph, false);

		    /* first find the first fundamental cycle back to the root
		     * do a normal dfs
		     */
		    std::cout << "first cycle " << std::endl;
		    if (!find_a_cycle_to_root(current_node, number_seen, seen_edge)) {
		    	assert(false); // graph is not triconnected
		    }

		    const node node_a = current_node;
		    /* When the first cycle is found we start a new DFS for every node on the path from a upwards to the root. When one
		     * of those searches encounters another backedge to the root we know that the starting node is the LCA.
		     *
		     * TODO we can also start labelling edges as belonging to a path
		     */

		    node inner_search_cur_node;
		    while(current_node != NULL) {
		    	std::cout << "walk up " << current_node->id() <<std::endl;

		    	inner_search_cur_node = current_node;
		    	if (find_a_cycle_to_root(inner_search_cur_node,number_seen, seen_edge)) {
		    		break;
		    	}

		    	current_node = parent[current_node];
		    }

		    const node lca = current_node;
		    const node node_b = inner_search_cur_node;

		    std::cout << "Root " << root->id() << std::endl;
		    std::cout << "A " << node_a->id() << std::endl;
		    std::cout << "B " << node_b->id() << std::endl;
		    std::cout << "LCA " << lca->id() << std::endl;

	}

	bool find_a_cycle_to_root(node& current_node, unsigned int& number_seen, edge_array<bool>& seen_edge) {
			stack<edge> first_cycle;

			{	edge e;
				forall_adj_edges(e,current_node) {
					if (seen_edge[e])
						continue;
					if (source(e)!=current_node)
						the_graph.rev_edge(e);


					parent[target(e)] = current_node;
					first_cycle.push(e);
				}
			}

			while(!first_cycle.empty()) {
				{
					const edge current_edge = first_cycle.pop();
					seen_edge[current_edge] = true;

					std::cout << source(current_edge)->id() << " " << target(current_edge)->id() << std::endl;

					current_node = source(current_edge);
					const node next_node = target(current_edge);

					current_node = next_node;
				}

				dfi[current_node] = number_seen++;

				{	edge next_edge;
					//push unvisited neighbours
					forall_adj_edges(next_edge,current_node) {

						//skip edge to parent
						if (seen_edge[next_edge])
							continue;



						if (source(next_edge)!=current_node)
							the_graph.rev_edge(next_edge);

						const node target_node = target(next_edge);
						// mark and skip fronds
						const unsigned int num = dfi[target_node];
						if (num != 0 && num<dfi[current_node]) {
							is_frond[next_edge] = true;
							if (num == 1) {
								//first cycle complete
								std::cout << "found cycle from " << current_node->id() << " " << target_node->id() << std::endl;
								return true;
							}
							continue; // just some other frond
						}

						parent[target_node] = current_node;

						first_cycle.push(next_edge);

					}
				}
			}

			return false;
	}
};

bool schmidt_is_triconnected(ugraph& g) {
	schmidt_triconnectivity t(g);
	return true;
}
