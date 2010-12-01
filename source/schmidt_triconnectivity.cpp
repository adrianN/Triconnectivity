#include "LEDA/graph/ugraph.h"
#include "LEDA/core/stack.h"
#include <fstream>
#include <iostream>

using namespace leda;
enum chain_type { one, two_a, two_b, three_a, three_b};

struct chain {
	chain* parent;
	node s;
	node t;
	chain_type type;

};





class schmidt_triconnectivity {


	edge_array<int> belongs_to_chain;
	edge_array<bool> is_frond;

	node_array<unsigned int> dfi;
	node_array<edge> parent;
	chain* chains;
	ugraph& the_graph;

public:
	schmidt_triconnectivity(ugraph& graph) :
		belongs_to_chain(graph,-1),
		is_frond(graph, false),
		dfi(graph,0),
		parent(graph, NULL),
		chains(new chain[graph.number_of_edges() - graph.number_of_nodes() +2]),
		the_graph(graph)
	{
		{	edge e;
			forall_edges(e,the_graph) {
				belongs_to_chain[e] = -1;
			}
		}
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
		    slist<edge> fronds;

		    chains[0].parent = NULL;
		    chains[1].parent = &chains[0];
		    chains[2].parent = &chains[0];

		    /* first find the first fundamental cycle back to the root
		     * do a normal dfs
		     */
		    std::cout << "first cycle " << std::endl;
		    if (!find_a_cycle_to_root(current_node, number_seen, seen_edge, false,  1, fronds)) {
		    	assert(false); // graph is not triconnected
		    }

		    const node node_a = current_node;
		    /* When the first cycle is found we start a new DFS for every node on the path from a upwards to the root. When one
		     * of those searches encounters another backedge to the root we know that the starting node is the LCA.
		     *
		     * TODO we can also start labelling edges as belonging to a path
		     */

		    node inner_search_cur_node;
			node lca = NULL;
		    node node_b = NULL;
		    while(current_node != NULL) { //todo switch from node to edge
		    	std::cout << "walk up " << current_node->id() <<std::endl;

		    	inner_search_cur_node = current_node;

		    	if (find_a_cycle_to_root(inner_search_cur_node,number_seen, seen_edge, true, 2, fronds) && lca == NULL && node_b == NULL) {
		    		std::cout << "found a cycle" << std::endl;
		    		lca = current_node;
		    		node_b = inner_search_cur_node;

		    	}
		    	if (parent[current_node]!=NULL) {
		    		std::cout << "the parent of " << current_node->id() << " is " << source(parent[current_node])->id() << " " << target(parent[current_node])->id() << std::endl;
					if (source(parent[current_node]) == current_node) {
						the_graph.rev_edge(parent[current_node]);
					}
					current_node = source(parent[current_node]);
		    	} else {
		    		current_node = NULL;
		    	}

		    }
		    assert(lca!=NULL);
		    assert(node_b!=NULL);
		    std::cout << "Root " << root->id() << std::endl;
		    std::cout << "LCA " << lca->id() << std::endl;

		    mark_path(root, lca, 0);
		    std::cout << "A " << node_a->id() << std::endl;

		    mark_path(root, node_a, 1);
		    std::cout << "B " << node_b->id() << std::endl;

		    mark_path(root, node_b, 2);


		    chain_decompose(root, fronds);
	}

	void mark_path(const node root, const node start, const unsigned int the_chain) {
		chains[the_chain].s = start;
		for(node current = start; current !=root; current = parent[current]!=NULL ? source(parent[current]) : NULL) {
			std::cout << "mark " << source(parent[current])->id() << " " << target(parent[current])->id() << std::endl;
			if (belongs_to_chain[parent[current]]>=0) {
				assert((unsigned int)belongs_to_chain[parent[current]] < the_chain);
				chains[the_chain].parent=&chains[belongs_to_chain[parent[current]]];
				break; //todo set father for chain
			}
			assert(dfi[chains[the_chain].s] >= dfi[current]);
			belongs_to_chain[parent[current]] = the_chain;
		}
	}

	void chain_decompose(const node root, slist<edge> fronds) {
		unsigned int chain_number = 3; //the first three chains are found during dfs
		assert(!fronds.empty());

		for(edge e=fronds.pop(); e!=NULL; e=(fronds.empty() ? NULL : fronds.pop())) {
			std::cout << dfi[source(e)] << " " << dfi[target(e)] << std::endl;
			assert(belongs_to_chain[e] <= 2);
			if (belongs_to_chain[e]<0) {
				belongs_to_chain[e] = chain_number;
				mark_path(root,source(e),chain_number);
				chain_number++;
			}
		}
	}

	bool find_a_cycle_to_root(node& start_node, unsigned int& number_seen, edge_array<bool>& seen_edge, bool continue_after_found, const unsigned int current_chain, slist<edge>& fronds) {
			stack<edge> first_cycle;
			node current_node = start_node;
			{	edge e;
				forall_adj_edges(e,current_node) {
					if (seen_edge[e])
						continue;

					assert(!seen_edge[e]);
					std::cout << "\tpush " << e << " " << source(e)->id() << " " << target(e)->id() << std::endl;

					first_cycle.push(e);
				}
			}

			bool found_cycle = false;
			while(!first_cycle.empty()) {
				{
					const edge current_edge = first_cycle.pop();

					//make sure fronds go from down to up and tree edges from up to down
					if (is_frond[current_edge]) {
						if (dfi[source(current_edge)] < dfi[current_node]) {
							the_graph.rev_edge(current_edge);
						}
						continue;
					} else {
						if (source(current_edge)!=current_node)
							the_graph.rev_edge(current_edge);
					}
					std::cout << source(current_edge)->id() << " " << target(current_edge)->id() << std::endl;

					if (seen_edge[current_edge]) {
						std::cout << "uh oh" << std::endl;
						std::cout.flush();
						assert(false);
					}
					seen_edge[current_edge] = true;


					current_node = source(current_edge);
					const node next_node = target(current_edge);
					assert(parent[next_node]==NULL);
					parent[next_node] = current_edge;
					std::cout << "\tparent " << next_node->id() << " " << source(current_edge)->id() << " " << target(current_edge)->id() << std::endl;

					current_node = next_node;
				}

				dfi[current_node] = number_seen++;

				{	edge next_edge;
					//push unvisited neighbours
					int deg = 0;
					std::cout << "\t deg " << the_graph.degree(current_node) << std::endl;
					forall_adj_edges(next_edge,current_node) {
						assert(++deg <= the_graph.degree(current_node));

						//skip edge to parent
						if (seen_edge[next_edge]) {
							std::cout << "\tskip " << next_edge << " " << source(next_edge)->id() << " " << target(next_edge)->id() <<  std::endl;
							continue;
						}

						const node target_node = opposite(next_edge,current_node);
						// mark fronds
						const unsigned int num = dfi[target_node];
						if (num != 0 && num<dfi[current_node]) {
							is_frond[next_edge] = true;
							fronds.append(next_edge);
							seen_edge[next_edge] = true;
							if (num == 1 && !found_cycle) { //backedge to the root
								//first cycle complete
								std::cout << "found cycle from " << current_node->id() << " " << target_node->id() << std::endl;
								belongs_to_chain[next_edge] = current_chain;
								start_node = current_node;
								found_cycle=true;
								if (!continue_after_found) {
									seen_edge[next_edge] = true;
									if (dfi[source(next_edge)] < dfi[current_node]) { //since we won't inspect this edge again we need to make sure it's the right way round here.
										std::cout << "reverse frond to root" << std::endl;
										the_graph.rev_edge(next_edge);
									}
									return found_cycle;
								}
							}

						}

						assert(!seen_edge[next_edge] || is_frond[next_edge]);
						std::cout << "\tpush " << next_edge << " " << current_node->id() << " " << target_node->id() << std::endl;

						first_cycle.push(next_edge);
					}
				}
			}

			return found_cycle;
	}



	   void to_dot(std::ostream& out) {
	        out << "digraph G {" << std:: endl;
	        {
	            node n;
	            forall_nodes(n,the_graph) {
	                unsigned int num = dfi[n];
	                int id = n->id();
	                /*
	                 * +--------------------------------+
	                 * | Parent             			|
	                 * +--------------------------------+
	                 * | LP1 | LP2 | highpoint | #desc. |
	                 * +--------------------------------+
	                 * | Node               			|
	                 * +--------------------------------+
	                 */
	                out << "\tnode";
	                out << id;
	                out << " [shape = record, label = \"{";
	                if (parent[n] == NULL) {
	                	out << "Father: " << 0;
	                } else {
	                	edge e = parent[n];

	                	out << "Father " << source(e);
	                }
	                out << "|Node " << num << ", ID " << id << "}\"]" << std::endl;

	            }
	        }
	        {
	        	/* Edges are labeled with 1 if they start a path (0 oth.) and the value they by which they get ordered (phi) */
	            edge e;
	            forall_edges(e,the_graph) {
	                node n = source(e);
	                node u = target(e);





	                out << "\tnode" << n->id() << " -> " << "node" << u->id();
	                out << " [ " << (is_frond[e] ? "constraint = false, color = \"red\", ":"") << "label=\"" << belongs_to_chain[e] <<  "\"]";
	                out << std::endl;
	            }

	        }


	        out << "}" << std::endl;
	    };

};

bool schmidt_is_triconnected(ugraph& g) {
	schmidt_triconnectivity t(g);
	std::fstream f("./blablubb.dot", std::ios::out);
	t.to_dot(f);
	return true;
}
