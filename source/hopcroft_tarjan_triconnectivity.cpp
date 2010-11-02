#include "triconnectivity.hpp"
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
		the_graph(g),
		number(the_graph,-1),
		number_descendants(the_graph,0),

		is_frond(the_graph,false),
		is_first_on_path(the_graph,0),
		lowpoint_one(new unsigned int[number_of_nodes+1]),
		lowpoint_two(new unsigned int[number_of_nodes+1]),
		flag(new bool[number_of_nodes+1]),
		father(new unsigned int[number_of_nodes+1]),
		highpoint(new int[number_of_nodes+1]),
		m_NODEAT(new node[number_of_nodes+1])

	{
		step_one();
		step_two();
		step_one();
	};

	~palm_tree() {
		delete[] lowpoint_one;
		delete[] lowpoint_two;
		delete[] father;
		delete[] highpoint;
		delete[] flag;
		delete[] m_NODEAT;
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

		unsigned int n = 0;
		dfs(the_graph.first_node(), 0, n);

		make_acceptable_adjacency();


	}

	unsigned int phi(const edge& e) const {
		const unsigned int v = number[source(e)];
		const unsigned int w = number[target(e)];
		if (is_frond[e]) return 3*w+1;
		if (lowpoint_two[w] < v) return 3*lowpoint_one[w];
		return 3*lowpoint_one[w]+1;
	}

	void make_acceptable_adjacency() {
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
	 *
	 * TODO check for biconnectedness
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
			if (node_w!=node_v)
				continue;

			int w = number[node_w];

			if (w < 0 || (w==0 && node_w!=the_graph.first_node())) { //latter part is fix for bug in LEDA
				dfs(node_w, v,n);

				w = number[node_w];

				if (lowpoint_one[w] < lowpoint_one[v]) {
					lowpoint_two[v] = min(lowpoint_one[v],lowpoint_two[w]);
					lowpoint_one[v] = lowpoint_one[w];
				} else if (lowpoint_one[w] == lowpoint_one[v]) {
					lowpoint_two[v] = min(lowpoint_two[v], lowpoint_two[w]);
				} else {
					lowpoint_two[v] = min(lowpoint_two[v],lowpoint_one[w]);
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

//	void TricComp::pathFinder(const Graph& G, node v)
//	{
//		m_NEWNUM[v] = m_numCount - m_ND[v] + 1;
//
//		ListConstIterator<edge> it;
//		for(it = m_A[v].begin(); it.valid(); ++it) {
//			edge e = *it;
//			node w = e->opposite(v);
//
//			if (m_newPath) {
//				m_newPath = false;
//				m_START[e] = true;
//			}
//
//			if (m_TYPE[e] == tree) {
//				pathFinder(G,w);
//				m_numCount--;
//
//			} else {
//				m_IN_HIGH[e] = m_HIGHPT[w].pushBack(m_NEWNUM[v]);
//				m_newPath = true;
//			}
//		}
//	}

	void pathfinder(node node_v, bool s, unsigned int m) {
		assert(m>=0 && m <=number_of_nodes);
		const int v = m - number_descendants[node_v]  +1 ;
		assert(v>=0);
		assert((unsigned int)v<number_of_nodes);
		number[node_v] = v;
		{ 	edge e;
			forall_adj_edges(e,node_v) {
				node node_w = target(e);
				if (node_w == node_v)
					continue;

				is_first_on_path[e] = s;
				if (s) s=false;

				if (!is_frond[e]) {
					pathfinder(node_w, s,m-1);
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

	bool pathsearch (node node_v, node &s1, node &s2)
	{

		unsigned int y;
		const unsigned int v = number[node_v];
		const node node_first_child = the_graph.adj_nodes(node_v).contents(the_graph.adj_nodes(node_v).first());
		unsigned int a, b;

		//List<edge> &Adj = m_A[v];
		unsigned int outv = the_graph.degree(node_v);

		edge e;
		forall_adj_edges(e,node_v) {
			//itNext = it.succ();
			node node_w = target(e);
			if (node_w == node_v) continue;
			const unsigned int w = number[node_w];

			if (!is_frond[e]) {
				if (is_first_on_path[e]) {
					y = 0;
					if (t_stack.top().second() > lowpoint_one[w]) {
						do {
							three_tuple<unsigned int, unsigned int, unsigned int> hab = t_stack.pop();
							y = max(y,hab.first());
							b = hab.second();
						} while (t_stack.top().second() > lowpoint_one[w]);
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(y,lowpoint_one[w],b));
					} else {
						t_stack.push(three_tuple<unsigned int, unsigned int, unsigned int>(w+number_descendants[node_w]-1,lowpoint_one[w], v));
					}
					//TSTACK_pushEOS();
				}

				if (!pathsearch(node_w,s1,s2))
					return false;

				while (v != 1 && ((t_stack.top().second() == v) ||
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
						s1 = m_NODEAT[a];
						s2 = m_NODEAT[b];
						return false;
					}
				}

				if (lowpoint_two[w] >= v && lowpoint_one[w] < v && (father[v] != (unsigned int)number[the_graph.first_node()] || outv >= 2))
				{
					s1 = m_NODEAT[lowpoint_one[w]];
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
					if (t_stack.top().second() > w) {
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
	ugraph& the_graph;
	node_array<int> number; //stores the step in which a node is reached during DFS
	node_array<unsigned int> number_descendants;


	edge_array<bool> is_frond; //true if edge is a frond.
	edge_array<bool> is_first_on_path;

	unsigned int* const lowpoint_one; //stores for each node the first lowpoint
	unsigned int* const lowpoint_two;
	bool* const flag; //becomes true when the entry in adj(v) corresponding to tree arc (u,v) is examined
	unsigned int* const father;
	int* const highpoint;
	node* const m_NODEAT;

	stack<three_tuple<unsigned int, unsigned int, unsigned int> > t_stack;
	stack<edge> e_stack;

};

void test(ugraph& g) {
	palm_tree p(g);
	std::fstream f("./the_palm_tree.dot", std::ios::trunc | std::ios::out);
	p.to_dot(f);
}
