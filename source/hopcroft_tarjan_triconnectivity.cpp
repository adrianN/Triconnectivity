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
		node_at(new node[number_of_nodes+1]),
		current_t_stack(new t_stack_type())

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

	typedef stack<three_tuple<unsigned int, unsigned int, unsigned int> > t_stack_type;
	auto_ptr<t_stack_type> current_t_stack;
	stack<auto_ptr<t_stack_type> > t_stacks;
	//stack<edge> e_stack;

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
			seen_edge_to_parent[i] = false;
			father[i] = 0xDEAD;
			lowpoint_one[i] = -1;
			lowpoint_two[i] = -1;
		}

		const int parent = -1;
		const int next_number = 1;
		dfs(the_graph.first_node(), parent, next_number,seen_edge_to_parent);

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


		if (is_frond[e]) return 3*w+1; // this goes in the middle. Remark: If the graph is biconnected, I don't think we need the +1 -- as lowpoint_one[w] < w --, but it doesn't hurt either.

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

		for(unsigned int i = 0; i<number_of_nodes; i++) {
			highpoint[i] = -1;
		}

		const bool at_start_of_path=true;
		const unsigned int nodes_to_be_seen = number_of_nodes;
		pathfinder(the_graph.first_node(), at_start_of_path, nodes_to_be_seen);

		the_graph.bucket_sort_nodes(number); //probably unnecessary

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
	 * Further each frond is marked as a frond in is_frond[e]
	 */
	void dfs(const node& node_v, const int parent, unsigned int next_number, bool* const seen_edge_to_parent) {
//		std::cout << "DFS " << node_v->id() << std::endl;

		number[node_v] = next_number++;

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
					dfs(node_w, v,next_number,seen_edge_to_parent);

					//during the dfs w gets assigned a number
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

				} else if ((unsigned int)w < v && ((w!=parent) || seen_edge_to_parent[v]) ) {
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
				if (w==parent) seen_edge_to_parent[v] = true;

				/* Here I check for biconnectedness. If a graph is biconnected, v -> w => lowpoint_one[w] < v, unless v=1. In that case lowpoint_one[w] = v = 1
				 * If lowpoint_one[w] >= v, we can cut v out of the graph and separate w (and its descendants) from the rest, since there is no path back beyond v
				 * This biconnectedness test is most likely not exhaustive. I'm waiting for a bugfix in LEDA TODO
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


	/* This function repeats the traversal of the pathfinder and looks for separation pairs in the graph. It stops after finding the first separation pair.
	 *
	 * (a,b) is a separation pair iff there are two other nodes x,y such that all paths between x and y cross a or b (or both); i.e. when a and b are removed, x and y are separated.
	 *
	 * There are two (or three) types of separation pairs (a,b), a<b
	 *
	 * 1. There are distinct vertices r != a,b and s!=a,b such that
	 * 		i) b->r,
	 * 		ii) lowpoint_one[r]=a,
	 * 		iii) lowpoint_two[r] >= b and
	 * 		iv) s is not a descendant of r
	 * 2. There is a vertex r!=b such that
	 * 		i) a->r-*>b,
	 * 		ii) b is a first descendant of r
	 * 		iii) a!=1
	 * 		iv) every frond x --> y with r<=x<b has a<=y
	 * 		v) and every frond x-->v with a<v<b and b -> w -*>x has lowpoint_one[w] >= a
	 * (3.) (a,b) is a multiple edge of G and G contains at least four edges
	 *
	 * =>
	 *
	 * Suppose (a,b) satisfies 1.; we need to show that there are at least two separation classes. The edge (b,r) is contained is some separation class, say E_1. Every tree arc starting in the subtree of
	 * r has the other end also there (paper says union {a,b}, but I don't see why those could be tree arcs, not fronds, as a,b<r). Since lowpoint_one[r]=a and lowpoint_two[r]>=b (ii, iii)
	 * every frond in the subtree either ends at a or stays in the subtree union b. Hence if a,b are removed the subtree under r is separated. There is another separation class E_2 because there
	 * still is s around, which isn't a descendant of r.
	 * Suppose (a,b) satisfies 2. Let S be the subtree of r - the subtree of b. All nodes in S are in one separation class, say E_1. Since b is a first descendant of r (ii), it gets finished first during
	 * the pathfinder dfs and has the highest number, hence S = {x| r<= x < b}. Let b_1,...,b_n be the sons of b. By the ordering in the acceptable adjacency structure, there is an index
	 * i_0 s.t. i < i_0 implies lowpoint_one[b_i] < a and i>=i_0 implies lowpoint_one[b_i] >= a. By (iv) every frond which starts in S ends in S union a, since a<=y and a->r there are no nodes
	 * y outside of S union a. (v) makes sure that fronds that reach into S from below b don't connect E_1 with nodes above a, that is, they start in S cup b cup D(b_i), i>=i_0. Thus E_1 contains
	 * at least all edges ending in S and at most all edges ending in S and D(b_i) i>=i_0 (if type (v) fronds exist). But since a is not the root of the tree, there must be some more nodes
	 * v < a, that are in a different separation class. That's why it is necessary to remove b too, as there could be fronds starting in the subtree of b that reach up to v<a.
	 *
	 * <=
	 *
	 *
	 * The t_stacks store possible candidates for type-2 separation pairs. The current_t_stack is the stack for the currently examined path. The stacks store triples (h,a,b); (a,b) is a
	 * possible type-2 separation pair and h is the highest number in the corresponding split component.
	 */
	bool pathsearch (node node_v, node &s1, node &s2)
	{
		const unsigned int v = number[node_v];

		unsigned int num_unexplored_adj_edges = the_graph.degree(node_v);

		edge e;
		forall_adj_edges(e,node_v) {
			num_unexplored_adj_edges--;

			if (is_frond[e]) continue;

			const node node_w = target(e);
			if (node_w == node_v) continue;
			const unsigned int w = number[node_w];

			if (is_first_on_path[e])
			current_t_stack	= update_stack(e,node_w, w, v);

			if (!pathsearch(node_w,s1,s2))
				return false;

			/* Check for type 2 pairs
			 * 2. There is a vertex r!=b such that
			 * i) a->r-*>b,
			 * ii) b is a first descendant of r
			 * iii) a!=1
			 * iv) every frond x --> y with r<=x<b has a<=y
			 * v) and every frond x-->v with a<v<b and b -> w -*>x has lowpoint_one[w] >= a
			 *
			 * With the current variable naming a = "v", r = "w".
			 *
			 * Things on the t_stack satisfy (i) and (ii) by the order of the search, (iv) by update_t_stack
			 */
			if (v!=1) {
				{	const node node_child_w = the_graph.adj_nodes(node_w).contents(the_graph.adj_nodes(node_w).first());
					// it should also be ok to check whether the first out-edge of w is a frond.
					// I don't know why we need to check number[node_child_w] > w. TODO
					const bool deg_two_case =  (the_graph.degree(node_w) == 2 && (unsigned int)number[node_child_w] > w);

					if (deg_two_case) {
						//separate node with degree two by its neighbours
						// v -- w -- child(w)
						s1 = node_v;
						s2 = node_child_w;
						return false;
					}
				}

				/* Iterate over all potential separation pairs with a = v */
				while (!current_t_stack->empty() && (current_t_stack->top().second() == v)) {
					const unsigned int 	b = current_t_stack->top().third();

					if (father[b] == v) { //this stack element is not a type two pair, as there are no nodes between a and b
						current_t_stack->pop();
					} else {
						s1 = node_v;
						s2 = node_at[b];
						return false;
					}
				}
			}

			/* Check for type one pairs
			 * 1. There are distinct vertices r != a,b and s!=a,b such that
			 * 		i) b->r,
			 * 		ii) lowpoint_one[r]=a,
			 * 		iii) lowpoint_two[r] >= b and
			 * 		iv) s is not a descendant of r
			 *
			 * With the current variable naming, r is "w", b is "v" and a is "lowpoint_one[w]". We need to make sure that there is a vertex s. For this it is sufficient
			 * to check that either father[v] is not the root (otherwise v > lowpoint_one[w] => lowpoint_one[w] = root = a) or v has other children than w.
			 *
			 * I'm not sure why it is correct to consider all neighbours of v including those connected by fronds.
			 */
			// i satisfied by choice of edge v -> w
			// if lowpoint_two[v] >= v and lowpoint_one[w] >= v removal of node v should disconnect the graph => not biconnected.
			if (/* iii */ lowpoint_two[w] >= v && /* a<b */lowpoint_one[w] < v && (father[v] != 1 || num_unexplored_adj_edges > 0))
			{
				s1 = node_at[lowpoint_one[w]]; // ii
				s2 = node_v;
				return false;
			}

			if (is_first_on_path[e]) {
//				if (t_stacks.empty())
//					current_t_stack->clear();
//				else
					current_t_stack = t_stacks.pop(); //auto_ptr should take care of itself here
			}

			/*
			 * This is another step for maintaining the t_stack.
			 */
			while (!current_t_stack->empty() && current_t_stack->top().third() != v && (unsigned int)highpoint[v] > current_t_stack->top().first()) {
				current_t_stack->pop();
			}


		}

		return true;
	}

	/* This function updates the t_stack. It puts all potential separation pairs on it (it's the only place where something gets pushed on the t_stack
	 * and removes those pairs that violate some conditions for type-2 pairs.
	 * the edge e starts a new path. Remember, a type-2 separation pair (a,b) isolates the part S between a and b from {the nodes above a, the subtree of b}.
	 *
	 * 2. There is a vertex r!=b such that
	 * i) a->r-*>b,
	 * ii) b is a first descendant of r
	 * iii) a!=1
	 * iv) every frond x --> y with x in S has a<=y
	 * v) and every frond x-->v with a<v<b and b -> w -*>x has lowpoint_one[w] >= a
	 *
	 * We make sure everything on the t_stack satisfies (i), (iv) and (v)
	 */
	auto_ptr<t_stack_type> update_stack(const edge e, const node node_w, const unsigned int w, const unsigned int v) {

		if (!is_frond[e]) {

			/* When we start a new path with the tree edge v -> w it's possible to violate (iv) for possible pairs (a,b) on the stack. So we check that lowpoint_one[w] <= a and pop,
			 * if necessary.
			 *
			 * If nothing was popped (lowpoint_one[w],v) is a possible separation pair (why? TODO). We push it on the stack, together with the highest node number in w's subtree.
			 * Otherwise, if some triple (h,a,b) was deleted last, we get a new possible sep. pair with (lowpoint_one[w],b), as this pair obviously satisfies (iv). We push that
			 */

			// in case nothing gets deleted, we want to push these values:

			const unsigned int highest_vertex_number = w + number_descendants[node_w] -1;
			three_tuple<unsigned int, unsigned int, unsigned int> hab(highest_vertex_number, lowpoint_one[w], v);

			// remove all tuples from the stack that violate (iv), remember the last one
			three_tuple<unsigned int, unsigned int, unsigned int> last_popped = hab;
			while(!current_t_stack->empty() && current_t_stack->top().second() > lowpoint_one[w]) {
				last_popped = current_t_stack->pop();
			}

			// (lowpoint_one[w],b) is a possible sep. pair
			hab.first() = max(highest_vertex_number,last_popped.first());
			hab.third() = last_popped.third();

			current_t_stack->push(hab);

			// Since we will continue with examining things on the new cycle, we save this t_stack for later, when we return.
			t_stacks.push(current_t_stack);
			return auto_ptr<t_stack_type>(new t_stack_type());

		} else {
			/* In the case of a frond, the new path consists of just one edge, so the top pair of the stack is ok, if w >= a (v). So we pop all triples with w<a.
			 *
			 * If nothing was deleted (w,v) is a possible type-2 pair. Since v is processed last in this split component it has the highest number and we push (v,w,v)
			 * Otherwise let (h,a,b) be the last triple that got deleted. Then (w,b) is a possible new pair. The highest number for that pair is the maximum that was deleted.
			 */

			unsigned int highest_vertex_number = v;
			three_tuple<unsigned int, unsigned int, unsigned int> hab(v,w,v);

			three_tuple<unsigned int, unsigned int, unsigned int> last_popped = hab;
			while(!current_t_stack->empty() && w < current_t_stack->top().second()) {
				last_popped = current_t_stack->pop();
				highest_vertex_number = max(highest_vertex_number, last_popped.first());
			}

			hab.first() = highest_vertex_number;
			hab.third() = last_popped.third();

			current_t_stack->push(hab);

			// Since this segment is done processing after this edge, we don't have to start a new t_stack

			return current_t_stack;
		}
	}
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
