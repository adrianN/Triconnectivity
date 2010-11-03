#include "triconnectivity.hpp"
#include "dfs.hpp"

#include "LEDA/graph/node_array.h"
#include "LEDA/graph/edge_array.h"
#include "LEDA/core/stack.h"
#include "LEDA/core/tuple.h"

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
		is_biconnected(true),
		the_graph(g),
		number(the_graph,-1),
		number_descendants(the_graph,0),

		is_frond(the_graph,false),
		is_first_on_path(the_graph,0),
		lowpoint_one(new unsigned int[number_of_nodes+1]),
		lowpoint_two(new unsigned int[number_of_nodes+1]),
		father(new unsigned int[number_of_nodes+1]),
		highpoint(new int[number_of_nodes+1]),
		node_at(new node[number_of_nodes+1])

	{
		step_one();
		step_two();
		step_one();
//		node s1,s2;
//		pathsearch(the_graph.first_node(),s1,s2);
	};

	~palm_tree() {
		delete[] lowpoint_one;
		delete[] lowpoint_two;
		delete[] father;
		delete[] highpoint;
		delete[] node_at;
	};

	bool is_triconnected(void) {
		node s1,s2;
		return pathsearch(the_graph.first_node(),s1,s2);
	}

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
				out << " [constraint = " << (is_frond[e] ? "false, color = \"red\"":"true") << ", label=\"" << is_first_on_path[e] << "\"]";
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

		bool* seen_edge_to_parent =new bool[number_of_nodes+1];
		for(unsigned int i=0; i<number_of_nodes; i++) {
			seen_edge_to_parent.get()[i] = false;
			father[i] = 0xDEAD;
			lowpoint_one[i] = -1;
			lowpoint_two[i] = -1;
		}

		dfs(the_graph.first_node(), -1, 1,seen_edge_to_parent);
		delete[] seen_edge_to_parent;
		make_acceptable_adjacency();


	}

	/* This function realizes the ordering of the edges. We want that
	 * - the edges vw are ordered ascending in order of lowpoint_one[w] (we want to create paths that end at the lowest possible point
	 * - for children w_1,...,w_n of v with lowpoint_one[w_i] = u (u<v) in the order given by the adjacency list. Then there exists an
	 *   i_0, such that lowpoint_two[w_i] < v for i<=i_0 and lowpoint_two[w_j] >= v for i_0 < j. If v --> u is in E then v --> u comes in Adj(v) between v -> w_i0 and v -> w_i0+1
	 *
	 * This way we make sure we always end up at the lowest possible point during the pathfinder, first traverse edges to nodes with a smaller lowpoint_one
	 */
	unsigned int phi(const edge& e) const {
		assert(number[source(e)] > 0);
		assert(number[target(e)] > 0);
		const unsigned int v = number[source(e)];
		const unsigned int w = number[target(e)];
		assert(w<=number_of_nodes);
		assert(v<=number_of_nodes);


		if (is_frond[e]) return 3*w+1; // this goes in the middle. Remark: If the graph is biconnected, I don't think we need the +1, but it doesn't hurt either.

		if (lowpoint_two[w] < v /* && !is_frond[e] implied */) return 3*lowpoint_one[w]; //lowpoint_one[w] < w in biconnected graphs, so this is smaller than the first number

		// if two children have the same lowpoint_one, order those with lowpoint_two > v at the end.
		return 3*lowpoint_one[w]+2;
	}

	/* order the edges according to phi, to make the adjacency structure acceptable */
	void make_acceptable_adjacency() {
		{	edge_array</*unsigned*/ int> phi_value(the_graph); //another bug in LEDA, unsigned ints are not integral.
			{	edge e;
				forall_edges(e,the_graph) {
					phi_value[e] = phi(e);
				}
			}

			the_graph.bucket_sort_edges(phi_value);
		}
		the_graph.bucket_sort_nodes(number); // this is probably unnecessary.
	}

	void step_two(void) {
		std::fstream f("./step_one.dot", std::ios::trunc | std::ios::out);
		to_dot(f);
		f.close();
		/*the smallest number
		 * Step two
		 */
		for(unsigned int i = 0; i<number_of_nodes; i++) {
			highpoint[i] = -1;
		}

		bool s=true;
		unsigned int m = number_of_nodes;
		pathfinder(the_graph.first_node(),s,m);

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

	/* This function performs a DFS on the graph and calculates a palm tree. It also calculates the following values
	 * - number[v], the order in which the vertices are seen during the search
	 * - lowpoint_one[v], the lowest node that can be reached by traversing >=0 tree edges and one frond
	 * - lowpoint_two[v], the _second_ lowest node that can be reached that way
	 * - number_descendants[v], the number of nodes in the subtree of v (including v)
	 * - father[v], the father of v in the palm tree
	 *
	 * Further each frond is marked as a frond in is_frond[e]the smallest number
	 */
	void dfs(const node& node_v, const int u, unsigned int n, bool* seen_edge_to_parent) {
//		std::cout << "DFS " << node_v->id() << std::endl;

		number[node_v] = n++;

		assert(number[node_v]>0);
		assert((unsigned int)number[node_v]<=number_of_nodes);
		const unsigned int v = number[node_v];

		// We initialize lowpoint one and two with v, in case no lower node can be reached
		lowpoint_one[v] = v;
		lowpoint_two[v] = v;

		// We know that v has itself in its subtree
		number_descendants[node_v] = 1;

		{	edge wv;
			forall_adj_edges(wv,node_v) {

				node node_w = target(wv);
				if (node_w==node_v) //we don't want to loop, LEDA seems to store each edge twice
					continue;

				int w = number[node_w];

				if (w < 0 || (w==0 && node_w!=the_graph.first_node())) { //latter part is fix for bug in LEDA
					//if w<0 we haven't seen w before, so we do a dfs on w
					dfs(node_w, v,n,seen_edge_to_parent);

					//during the dfs w get's assigned a number
					w = number[node_w];
					assert(w>0);
					assert((unsigned int)w == v+1);

					father[w] = v;

					//now we update lowpoints. If we can reach a lower node from w, we can also reach a lower node from v
					if (lowpoint_one[w] < lowpoint_one[v]) {
						// the lowest reachable node from w is lower than from v. We save that as the new lowpoint_one[v]. We also check for the second lowest node[w].
						// If it is lower than lowpoint_one[v] we use it as lowpoint_two[v]
						lowpoint_two[v] = min(lowpoint_one[v],lowpoint_two[w]);
						lowpoint_one[v] = lowpoint_one[w];
					} else if (lowpoint_one[w] == lowpoint_one[v]) {
						// if the lowpoints are the same we need to compare the second lowest reachable node
						lowpoint_two[v] = min(lowpoint_two[v], lowpoint_two[w]);
					} else {
						// lowpoint_one[w] isn't lower than lowpoint_one[v], but it could be lower than lowpoint_two[v]
						lowpoint_two[v] = min(lowpoint_two[v],lowpoint_one[w]);
					}

					// add all the nodes from w's subtree to v's descendants
					number_descendants[node_v] += number_descendants[node_w];

				} else if ((unsigned int)w < v && ((w!=u) || seen_edge_to_parent[v]) ) {
					/* In this case we have already seen w. We need to check whether the back-edge qualifies as a frond.
					 * An edge v --> w is a frond if w -*> v is a path in the tree (and v -- w -* v is a circle in the graph).
					 * So edges back to the father node u are not fronds, unless there are more edges u--v
					 *
					 * (it is possible to see edges with w>v
					 * 	v - x     1 - 2
					 *  |  /      |  /
					 *  | /       | /
					 *   w         3
					 *  )
					 */

					is_frond[wv] = true;

					/*
					 * A frond allows us to go back up in the graph, so we need to update the lowpoints
					 */
					if ((unsigned int)w<lowpoint_one[v]) {
						lowpoint_two[v] = lowpoint_one[v];
						lowpoint_one[v] = w;
					} else if ((unsigned int)w>lowpoint_one[v]) {
						lowpoint_two[v] = min(lowpoint_two[v],(unsigned int)w);
					}

				}

				// If we see another edge to u, it's a frond, as there is a path .
				if (w==u) seen_edge_to_parent[v] = true;

				/* Here I check for biconnectedness. If a graph is biconnected, v -> w => lowpoint_one[w] < v, unless v=1. In that case lowpoint_one[w] = v = 1
				 * If lowpoint_one[w] >= v, we can cut v out of the graph and separate w (and its descendants) from the rest, since there is no path back beyond v
				 */
				if (v!=1) {
					is_biconnected &= lowpoint_one[w] < v;
				} else {
					is_biconnected &= lowpoint_one[w] == v;
				}
			}
		}
	};

	/**
	 * This function does a dfs on an acceptable adjacency structure and finds special edge disjoint paths in the palm tree.
	 * A path as found by this function always ends at a frond and, because of the ordering in the adjecency structure,
	 * always ends at the lowest possible vertex.
	 *
	 * It also numbers the nodes in the inverse order of traversal, according to the ordering in the acceptable adjacency structure, and computes
	 * - node_at[v], a reverse map from node numbers to LEDA nodes
	 * - highpoint[v], for nodes at which fronds end, the node with the lowest (according to the new numbering) number at the source of a frond.
	 */
	void pathfinder(node node_v, bool at_path_start, unsigned int nodes_to_be_seen) {
		assert(nodes_to_be_seen <=number_of_nodes);
		std::cout << "Pathfinder " << nodes_to_be_seen << " " << number_descendants[node_v] << std::endl;

		/* because we want the numbers in reverse, but need them before we come back (for the highpoints) we count downwards.
		 * so we decrease the number each time we discover a vertex. Since we have number_descendants[v]-1 vertices to
		 * see when we first visit v, we will decrease m by that amount.
		 */
		assert(nodes_to_be_seen-(number_descendants[node_v]-1) > 0);
		const unsigned int v = nodes_to_be_seen - (number_descendants[node_v] - 1);
		assert(v<=number_of_nodes);
		number[node_v] = v;
		node_at[v] = node_v;

		{ 	edge vw;
			forall_adj_edges(vw,node_v) {
				/* In this loop we check whether the edge vw is the first on one of the paths */

				const node node_w = target(vw);
				if (node_w == node_v) // again this is to handle LEDA's peculiarities
					continue;

				is_first_on_path[vw] = at_path_start;
				if (at_path_start) at_path_start=false;

				if (!is_frond[vw]) {
					pathfinder(node_w, at_path_start, nodes_to_be_seen-1);
				} else {
					/* if we encounter a frond, it is the last edge on the path, so set s to true to mark the next edge as a start of a path.
					 * we also update the highpoint. We don't need to check whether it is smaller than before because of the ordering in which we traverse the graph
					 */
					const unsigned int w = number[node_w];
					if (highpoint[v] == -1) {
						highpoint[w] = v;
					}
					at_path_start = true;
				}
			}
		}
	}


	bool pathsearch (node node_v, node &s1, node &s2)
	{
		unsigned int y;
		const unsigned int v = number[node_v];
		const node node_first_child = the_graph.adj_nodes(node_v).contents(the_graph.adj_nodes(node_v).first());
		unsigned int a, b;

		unsigned int outv = the_graph.degree(node_v);

		edge e;
		forall_adj_edges(e,node_v) {
			node node_w = target(e);
			if (node_w == node_v) continue;
			const unsigned int w = number[node_w];

			if (!is_frond[e]) {
				if (is_first_on_path[e]) {
					y = 0;
					if (!t_stack.empty() && t_stack.top().second() > lowpoint_one[w]) {
						do {
							three_tuple<unsigned int, unsigned int, unsigned int> hab = t_stack.pop();
							y = max(y,hab.first());
							b = hab.second();
						} while (t_stack.top().second() > lowpoint_one[w]);
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(y,lowpoint_one[w],b));
					} else {
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(w+number_descendants[node_w]-1,lowpoint_one[w], v));
					}
				}

				if (!pathsearch(node_w,s1,s2))
					return false;

				while (v != 1 && ((!t_stack.empty() && (t_stack.top().second() == v)) ||
					(the_graph.degree(node_w) == 2 && (unsigned int)number[node_first_child] > w)))
				{
					a = t_stack.top().second();
					b = t_stack.top().third();

					if (a == v && father[b] == a) {
						t_stack.pop();
					} else if (the_graph.degree(node_w) == 2 && (unsigned int)number[node_first_child] > w) //separate node with degree two by its neighbours
					{
						s1 = node_v;
						s2 = node_first_child;
						return false;

					} else {
						s1 = node_at[a];
						s2 = node_at[b];
						return false;
					}
				}

				if (lowpoint_two[w] >= v && lowpoint_one[w] < v && (father[v] != (unsigned int)number[the_graph.first_node()] || outv >= 2))
				{
					s1 = node_at[lowpoint_one[w]];
					s2 = node_v;
					return false;
				}

				if (is_frond[e]) {
					t_stack.clear();
				}

				while (!t_stack.empty() &&
					t_stack.top().third() != v && (unsigned int)highpoint[v] > t_stack.top().first()) {
					t_stack.pop();
				}

				outv--;

			} else { // frond arc
				if (is_first_on_path[e]) {
					y = 0;
					if (!t_stack.empty() && t_stack.top().second() > w) {
						do {
							y = max(y,t_stack.top().first());
							b = t_stack.pop().third();
						} while (t_stack.top().second() > w);
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(y,w,b));
					} else {
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(v,w,v));
					}
				}
			}
		}

		return true;
	}

	const unsigned int number_of_nodes;
	bool is_biconnected;
	ugraph& the_graph;

	node_array<int> number; //stores the step in which a node is reached during DFS
	node_array<unsigned int> number_descendants;


	edge_array<bool> is_frond; //true if edge is a frond.
	edge_array<bool> is_first_on_path;

	unsigned int* const lowpoint_one; //stores for each node the first lowpoint
	unsigned int* const lowpoint_two;
	unsigned int* const father;
	int* const highpoint;
	node* const node_at;

	stack<three_tuple<unsigned int, unsigned int, unsigned int> > t_stack;
	stack<edge> e_stack;

};

bool hopcroft_tarjan_is_triconnected(const ugraph& g) {
	if (!is_connected(g))
		return false;

	ugraph H;
	CopyGraph(H,g);
	palm_tree p(H);
	return p.is_triconnected();
}

void test(ugraph& g) {
	palm_tree p(g);
	std::fstream f("./the_palm_tree.dot", std::ios::trunc | std::ios::out);
	p.to_dot(f);
}
