#include "triconnectivity.hpp"
#include "LEDA/graph/node_array.h"
#include "LEDA/graph/edge_array.h"
#include <fstream>
#include <ostream>
#include <memory>

using std::auto_ptr;

using namespace leda;


/*
 * Finding separation pairs as described in "Dividing a graph into triconnected components"
 */

class palm_tree {
public:

	palm_tree(ugraph& g) :
		number_of_nodes((assert(g.number_of_nodes()>=0), g.number_of_nodes())),
		the_graph(g),
		number(the_graph,-1),
		number_descendants(the_graph,0),

		is_frond(the_graph,false),
		belongs_to_path(the_graph,0),
		lowpoint_one(new unsigned int[number_of_nodes]),
		lowpoint_two(new unsigned int[number_of_nodes]),
		flag(new bool[number_of_nodes]),
		father(new unsigned int[number_of_nodes]),
		highpoint(new int[number_of_nodes])

	{
		step_one();
		//step_two();
		//step_one();
	};

	~palm_tree() {
		delete[] lowpoint_one;
		delete[] lowpoint_two;
		delete[] father;
		delete[] highpoint;
		delete[] flag;
	};

	void to_dot(std::ostream& out) {
		out << "digraph G {" << std:: endl;
		{
			node n;
			forall_nodes(n,the_graph) {
				int num = number[n];
				/*
				 * +--------------------+
				 * | Parent             |
				 * +--------------------+
				 * | LP1 | LP2 | #desc. |
				 * +--------------------+
				 * | Node               |
				 * +--------------------+
				 */
				out << "\tnode";
				out << num;
				out << " [shape = record, label = \"{";
				out << "Father: " << father[num];
				out << " | {" << lowpoint_one[num];
				out << " | " << lowpoint_two[num];
				out << " | " << number_descendants[n];
				out << " } |Node " << num << ", ID " << n->id() << "}\"]" << std::endl;
			}
		}
		{
			edge e;
			forall_edges(e,the_graph) {
				unsigned int v,w;
				const node n = source(e);
				const node u = target(e);
				if (number[n] < number[u]) {
					v = number[n];
					w = number[u];
				} else {
					v = number[u];
					w = number[n];
				}
				if (is_frond[e]) {
					swap(v,w);
				}
				out << "\tnode" << v << " -> " << "node" << w;
				//assert(is_frond[e] || father[w] == v);
				out << " [constraint = " << (is_frond[e] ? "false, color = \"red\"":"true") << ", label=\"" << belongs_to_path[e] << "\"]";
				out << std::endl;
			}
		}


		out << "}" << std::endl;
	};

private:

	void step_one(void) {
		/*
		 * Step one, compute information and reorder edges
		 */

		//we depend on flag for conditionals
		for(unsigned int i=0; i<number_of_nodes; i++) {
			flag[i] = false;
			father[i] = 0xDEAD;
			lowpoint_one[i] = -1;
			lowpoint_two[i] = -1;
		}
		{	node n;
			forall_nodes(n,the_graph) {
				number[n] = -1;
				number_descendants[n] = -1;
			}
		}
		{	edge e;
			forall_edges(e,the_graph) {
				is_frond[e] = false;
			}
		}

		std::cout << "DFS from id " << the_graph.first_node()->id() << std::endl;
		unsigned int n = 0;
		dfs(the_graph.first_node(), -1, n);


		{	edge_array</*unsigned*/ int> phi_value(the_graph); //another bug in LEDA, unsigned ints are not integral.
			{	edge e;
				forall_edges(e,the_graph) {
					phi_value[e] = phi(e);
				}
			}

			the_graph.bucket_sort_edges(phi_value);
		}
		the_graph.bucket_sort_nodes(number);


	}

	void step_two(void) {
		std::fstream f("./step_one.dot", std::ios::trunc | std::ios::out);
		to_dot(f);
		f.close();
		/*
		 * Step two
		 */
		for(unsigned int i = 0; i<number_of_nodes; i++) {
			highpoint[i] = -1;
		}

		unsigned int current_path=0;
		bool s=true;
		unsigned int m = number_of_nodes;
		pathfinder(the_graph.first_node(),current_path,s,m);

		the_graph.bucket_sort_nodes(number);
		f.open("./step_2.dot", std::ios::trunc | std::ios::out);
		to_dot(f);
	}

	/* Step 1 in Hopcroft Tarjan algorithm. Compute the graph information.
	 * Things changed:
	 * - flag switched
	 * - instead of looking at adjacent nodes I traverse adjacent edges, to be able to mark edges as fronds
	 *
	 * The second argument is an int because parents already have a number.
	 */
	void dfs(const node& node_v, const int u, unsigned int& n) {

		number[node_v] = n++;

		assert(number[node_v]>=0);
		assert((unsigned int)number[node_v]<number_of_nodes);
		const unsigned int v = number[node_v];

		//addition a
		lowpoint_one[v] = v;
		lowpoint_two[v] = v;
		number_descendants[node_v] = 1; //why a one here? we don't know yet if there are any descendants. Maybe this also counts v as its zeroth descendant

		edge wv;
		forall_adj_edges(wv,node_v) {
			node node_w = target(wv);
			if (node_w == node_v) node_w = source(wv);
			assert(node_v != node_w);

			int w = number[node_w];

			if (w < 0 || (w==0 && node_w!=the_graph.first_node())) { //latter part is fix for bug in LEDA
				dfs(node_w, v,n);

				w = number[node_w];

				if (lowpoint_one[w] < lowpoint_one[v]) {
					lowpoint_two[v] = min(lowpoint_one[v],lowpoint_two[w]); //(lowpoint_one[v] < lowpoint_two[w]) ? lowpoint_one[v] : lowpoint_two[w];
					lowpoint_one[v] = lowpoint_one[w];
				} else if (lowpoint_one[w] == lowpoint_one[v]) {
					lowpoint_two[v] = min(lowpoint_two[v], lowpoint_two[w]); //(lowpoint_two[v] < lowpoint_two[w]) ? lowpoint_two[v] : lowpoint_two[v];
				} else {
					lowpoint_two[v] = min(lowpoint_two[v],lowpoint_one[w]); //(lowpoint_two[v] < lowpoint_one[w]) ? lowpoint_two[v] : lowpoint_one[w];
				}
				number_descendants[node_v] += number_descendants[node_w];
				father[w] = v;

			} else if ((unsigned int)w < v && ((w!=u) || flag[v]) ) {
				is_frond[wv] = true;

				if ((unsigned int)w<lowpoint_one[v]) {
					lowpoint_two[v] = lowpoint_one[v];
					lowpoint_one[v] = w;
				} else if ((unsigned int)w>lowpoint_one[v]) {
					lowpoint_two[v] = min(lowpoint_two[v],(unsigned int)w); //(lowpoint_two[v] < (unsigned int)w) ? lowpoint_two[v] : w;
				}
			}

			if (w==u) flag[v] = true;
		}
	};

	void pathfinder(node node_v, unsigned int& current_path, bool& s, unsigned int& m) {
		assert(m>=0 && m <=number_of_nodes);
		int v = m - number_descendants[node_v] /* +1 */; //unsure if +1, we start numbering with zero
		assert(v>=0);
		assert((unsigned int)v<number_of_nodes);
		number[node_v] = v;
		{
			{ 	edge e;
				forall_adj_edges(e,node_v) {
					if (belongs_to_path[e]!=0) continue;
					node node_w = target(e);
					if (node_w==node_v) node_w = source(e);

					if (s) {
						s = false;
						current_path++;
					}

					belongs_to_path[e] = current_path;

					if (!is_frond[e]) {
						pathfinder(node_w, current_path, s,m);
						m--;
					} else {
						unsigned int w = number[node_w];
						if (highpoint[v] == -1) {
							highpoint[w] = v;
						}
						s = true;
					}
				}
			}
		}

	}

	unsigned int phi(const edge& e) const {
		const unsigned int v = number[source(e)];
		const unsigned int w = number[target(e)];
		if (is_frond[e]) return 2*w+1;
		if (lowpoint_two[w] < v) return 2*lowpoint_one[w];
		return 2*lowpoint_one[w]+1;
	}

	const unsigned int number_of_nodes;
	ugraph& the_graph;
	node_array<int> number; //stores the step in which a node is reached during DFS
	node_array<unsigned int> number_descendants;


	edge_array<bool> is_frond; //true if edge is a frond.
	edge_array<unsigned int> belongs_to_path;

	unsigned int* const lowpoint_one; //stores for each node the first lowpoint
	unsigned int* const lowpoint_two;
	bool* const flag; //becomes true when the entry in adj(v) corresponding to tree arc (u,v) is examined
	unsigned int* const father;
	int* const highpoint;

};

void test(ugraph& g) {
	palm_tree p(g);
	std::fstream f("./the_palm_tree2.dot", std::ios::trunc | std::ios::out);
	p.to_dot(f);
}
