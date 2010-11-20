#include "triconnectivity.hpp"
#include "dfs.hpp"
#include "utilities.hpp"

#include "LEDA/graph/node_array.h"
#include "LEDA/graph/edge_array.h"
#include "LEDA/core/stack.h"
#include "LEDA/core/tuple.h"

#include <fstream>
#include <ostream>
#include <memory>

using std::auto_ptr;

using namespace leda;

#define DOTS
//#define DFSCOUT
//#define PHI_COUT
//#define PATHFINDER_COUT
#define PATHFINDER_PATH_COUT
#define PATHSEARCH_COUT
#define STACK_COUT
#define PATHSEARCH_COUT_SEPPAIR
#define RECORD_NODES

/*
 * Finding separation pairs as described in "Dividing a graph into triconnected components"
 *
 * It finds a separation pair of G if G is not triconnected.
 *
 * Let (a,b) be a pair of vertices in the biconnected multigraph G. Suppose the edges of G
 * are divided into equivalence classes E1..En such that two edges are in the same class
 * iff there is a path that contains both and doesn't contain a or b, except as an endpoint.
 * If there are at least two Ei, then {a,b} is a separation pair of G, unless
 * 	i) 	there are exactly two separation classes and one class consists of a single edge (this
 * 		only happens if the edge between (a,b) is a class by itself or G is not biconnected
 * 		and a or b is an articulation point)
 * 	ii) there are exactly three classes, each consisting of a single edge. (TODO think of an
 * 		example for this one that isn't x - a - b - y)
 *
 * Note that the definition of separation pair makes any pair of vertices a - b with a multiple
 * edge between them a separation pair, if there is at least one more edge in the graph.
 *
 * The algorithm is based on the idea of finding cycles in biconnected graphs and the removing
 * them to obtain a number of segments
 *
 * Let G be a biconnected multigraph and let c be a cycle in G. Let S1...Sn be the connected
 * components of G-c. S1 .. Sn and c partition the edges of G. Let {a,b} be a separation pair of
 * G such that (a,b) is not a multiple edge. Then the following things hold
 *
 * i) Either a,b is contained on c, or (a,b) is contained in _one_ segment
 *
 * Let, for the sake of contradiction, (a,b) be in segments S1, S2. Since G is biconnected
 * 		S1 and S2 are connected by two paths, both on c. For (a,b) to be a separation pair it must
 * 		induce at least two equivalence classes on the edges. All the edges in S1 are in one
 * 		equivalence class, say E1, since S1 is biconnected and a can't be an articulation point. The same
 * 		holds for S2. Clearly all other segments are also in E1 or E2, as they are connected
 * 		to S1, respectively S2, via >=2 paths, at least one of which doesn't contain a or b. But since
 * 		S1 and S2 are themselves connected to c by at least two edges each, E_1 = E_2 and (a,b)
 * 		is not a separation pair.
 *
 * ii) Suppose a and b both lie on c. Let p,q be the two paths on c that connect a with b. Then either
 *
 * 	a) 	Some segment Si with at least two edges has only a and b in common with c, and some
 * 		vertex v does not lie in Si (type 1 pair)
 * 	b) 	no segment contains a vertex v != a,b in p and a vertex w != a,b on q and p and q each contain a
 * 		vertex besides a and b (type 2 pair)
 *
 *		Suppose neither a) nor b) hold.
 *		If a) does not hold every segment has more nodes in common with c or there is no vertex that
 *		doesn't lie in Si. If there is no vertex outside of Si, (a,b) is not a separation pair
 *
 *
 * The idea of the algorithm is to find cycles c with their corresponding segments and check
 * (recursively) for the above conditions. Since both enumerating all cycles and checking all
 * pairs on the cycles in too expensive, something clever has to be done.
 *
 * We can use a DFS to find the cycles, but to simulate the recursive nature of the algorithm, we
 * need that
 *  a) the first DFS path is a cycle
 *  b) any consecutive path doesn't contain edges that lie on previous cycles and can be
 *  	completed to a cycle by adding part of exactly one old path (this way the cycle lies within
 *  	one segment.
 *
 * Since the DFS order depends on the order in the adjacency lists, we run a preliminary computation
 * that reorders the lists and builds an acceptable adjacency structure. This is done in step_one()
 *
 * We then proceed to identify the paths in the graph (step_two()) and finally, in another DFS, check
 * for separation pairs.
 *
 */

class palm_tree {
public:

    palm_tree(ugraph& g) :
        number_of_nodes((assert(g.number_of_nodes()>=0), g.number_of_nodes())),
        is_biconnected(true),
        the_graph(g),
        number(the_graph,0),
        number_descendants(the_graph,0),
        highpoint(the_graph,0),

        is_frond(the_graph,false),
        is_first_on_path(the_graph,false),
        lowpoint_one(new unsigned int[number_of_nodes+1]),
        lowpoint_two(new unsigned int[number_of_nodes+1]),
        father(new unsigned int[number_of_nodes+1]),
        node_at(new node[number_of_nodes+1]),
        current_t_stack(new t_stack_type())

    {
#ifdef DOTS
    	num_call=0;
        std::cout << "Dot output" << std::endl;
        {	std::fstream initial("./graph_initial.dot",std::ios::out);
			simple_to_dot(the_graph,initial);
			initial.close();
        }
#endif
#ifdef DFSCOUT
        std::cout << "Output during DFS" << std::endl;
#endif
#ifdef PATHFINDER_COUT
        std::cout << "Output during Pathfinder" << std::endl;
#endif
#ifdef STACK_COUT
        std::cout << "Output during Stackupdates" << std::endl;
#endif
#ifdef PATHSEARCH_COUT
        std::cout << "Extended Output during Pathsearch" << std::endl;
#endif
#ifdef PATHSEARCH_COUT_SEPPAIR
        std::cout << "Output during Pathsearch" << std::endl;
#endif


        step_one();
        step_two();


    };

    ~palm_tree() {
        delete[] lowpoint_one;
        delete[] lowpoint_two;
        delete[] father;
        delete[] node_at;
        delete current_t_stack;
        while(!t_stacks.empty()) {
        	delete t_stacks.pop();
        }
    };

    /* returns true if the graph is triconnected. Otherwise s1,s2 contain a
     * separation pair. If g is not even biconnected s1==s2 is an articulation point */
    bool is_triconnected(node& s1, node& s2) {

#ifdef DOTS
        std::fstream f("./triconnectivity.dot", std::ios::trunc | std::ios::out);
        to_dot(f);
        f.close();
#endif

    	if (!is_biconnected) {
#ifdef PATHSEARCH_COUT_SEPPAIR
    		std::cout << "Graph is not biconnected, articulation point " << articulation_point->id() << std::endl;
#endif
    		s1 = articulation_point;
    		s2 = articulation_point;
    		return false;
    	}
        return pathsearch(the_graph.first_node(),s1,s2);
    }

    void to_dot(std::ostream& out) {
        out << "digraph G {" << std:: endl;
        {
            node n;
            forall_nodes(n,the_graph) {
                int num = number[n];
                int id = n->id();
#ifdef RECORD_NODES
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
                out << "Father: " << father[num];
                out << " | {" << lowpoint_one[num];
                out << " | " << lowpoint_two[num];
                out << " | " << highpoint[n];
                out << " | " << number_descendants[n];
                out << " } |Node " << num << ", ID " << id << "}\"]" << std::endl;
#else
                out << "\tnode";
                out << id;
                out << " [label = \"";
                out << "Father: " << father[num] << "\\n";
                out << "LP1 " << lowpoint_one[num] << "\\n";
                out << "LP2 " << lowpoint_two[num] << "\\n";
                out << "HP " << highpoint[num] << "\\n";
                out << "ND " << number_descendants[n] << "\\n";
                out << "Node " << num << ", ID " << id <<  "\"]" << std::endl;
#endif
            }
        }
        {
            edge e;
            unsigned int num_of_fronds = 0, paths=0;
            forall_edges(e,the_graph) {
                node n = source(e);
                node u = target(e);

                if (number[n] > number[u]) {
                    swap(n,u);
                }
                if (is_frond[e]) {
                	num_of_fronds++;
                    swap(n,u);
                }
                if (is_first_on_path[e]) {
                	paths++;
                }
                out << "\tnode" << n->id() << " -> " << "node" << u->id();
                out << " [ " << (is_frond[e] ? "constraint = false, color = \"red\", ":"") << "label=\"" << is_first_on_path[e] << " " << phi(e) <<  "\"]";
                out << std::endl;
            }

            out << "label = \"Fronds: " << num_of_fronds << ", Paths " << paths << ", Nodes " << number_of_nodes << "\"" << std::endl;
        }


        out << "}" << std::endl;
    };

private:

    const unsigned int number_of_nodes;
    bool is_biconnected;
    node articulation_point;
    ugraph& the_graph;

    node_array<unsigned int> number; //stores the step in which a node is reached during DFS
    node_array<unsigned int> number_descendants;
    node_array<unsigned int> highpoint;



    edge_array<bool> is_frond; //true if edge is a frond.
    edge_array<bool> is_first_on_path;

    unsigned int* lowpoint_one; //stores for each node the first lowpoint
    unsigned int* lowpoint_two;
    unsigned int* father;
    node* const node_at;

    typedef stack<three_tuple<unsigned int, unsigned int, unsigned int> > t_stack_type;
    t_stack_type* current_t_stack;
    stack<t_stack_type* > t_stacks;
    //stack<edge> e_stack;

#ifdef DOTS
	int num_call;
#endif

    void step_one(void) {
#ifdef DOTS
    	num_call++;
#endif
        /*
         * Step one, compute information and reorder edges
         */

        {   node n;
            forall_nodes(n,the_graph) {
                number[n] = 0;
                number_descendants[n] = 0;
            }
        }
        {   edge e;
            forall_edges(e,the_graph) {
                is_frond[e] = false;
            }
        }

        bool* seen_edge_to_parent =new bool[number_of_nodes+1];
        for(unsigned int i=0; i<number_of_nodes+1; i++) {
            seen_edge_to_parent[i] = false;
            father[i] = 0;
            lowpoint_one[i] = 0; // but it is used that this is <0
            lowpoint_two[i] = 0;
        }

        const unsigned int parent = 0;
        const unsigned int next_number = 1;
        dfs_one(the_graph.first_node(), parent, next_number,seen_edge_to_parent);

        delete[] seen_edge_to_parent;
#ifdef DOTS
        {
		std::ostringstream s;
		s << "./dfs" << num_call << ".dot";
		std::fstream f(s.str().c_str(), std::ios::out);
        std::cout << " dfs dot " << std::endl;
        to_dot(f);
        f.close();
        }
#endif

        make_acceptable_adjacency();

#ifdef DOTS
        {
		std::ostringstream s;
		s << "./acceptable_adjacency" << num_call << ".dot";
		std::fstream f(s.str().c_str(), std::ios::out);
        to_dot(f);
        f.close();
        }
#endif

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
        unsigned int v = number[source(e)];
        unsigned int w = number[target(e)];

        if (is_frond[e]) { //fronds go from big to small, other nodes from small to big.
        	if (v<w) swap(v,w);
        } else {
        	if (v>w) swap(v,w);
        }


        assert(w<=number_of_nodes);
        assert(v<=number_of_nodes);

#ifdef PHI_COUT
        unsigned int phi = -1;
        unsigned int type = -1;

        if (is_frond[e]) {phi = 3*w+1; type = 1;}
        else if (lowpoint_two[w] < v) {phi = 3*lowpoint_one[w]; type = 2;}//lowpoint_one[w] < w in biconnected graphs, so this is smaller than the first number
        else {phi = 3*lowpoint_one[w]+2; type=3;}

        std::cout << "Edge-phi " << v << " " << w << " " << phi << " " << type << std::endl;
#endif

        if (is_frond[e]) return 3*w+1; // this goes in the middle. Remark: If the graph is biconnected, I don't think we need the +1 -- as lowpoint_one[w] < w --, but it doesn't hurt either.

        if (lowpoint_two[w] < v) return 3*lowpoint_one[w]; //lowpoint_one[w] < w in biconnected graphs, so this is smaller than the first number

        // if two children have the same lowpoint_one, order those with lowpoint_two > v at the end.
        return 3*lowpoint_one[w]+2;
    }

    /* order the edges according to phi, to make the adjacency structure acceptable */
    void make_acceptable_adjacency() {
        {   edge_array</*unsigned*/ int> phi_value(the_graph); //another bug in LEDA, unsigned ints are not integral.
            {   edge e;
                forall_edges(e,the_graph) {
                    phi_value[e] = phi(e);
                }
            }

            the_graph.bucket_sort_edges(phi_value);
        }
        //the_graph.bucket_sort_nodes(number); // this is probably unnecessary.
    }

    void step_two(void) {
    	unsigned int* old_new_map;
    	{	node_array<unsigned int> new_number(the_graph,0);

			{ node n; //BUGFIX for LEDA
				forall_nodes(n,the_graph) {
					highpoint[n] = 0;
					new_number[n] = 0;
				}
			}
	//        for(unsigned int i = 0; i<number_of_nodes+1; i++) {
	//            highpoint[i] = 0; //this is used
	//        }
			{	edge e; //BUGFIX for LEDA
				forall_edges(e,the_graph) {
					is_first_on_path[e] = false;
				}
			}

			bool at_start_of_path=true;
			const unsigned int nodes_to_be_seen = number_of_nodes;

			pathfinder(the_graph.first_node(), at_start_of_path, nodes_to_be_seen,new_number);

			old_new_map = new unsigned int[number_of_nodes+1];
			old_new_map[0] = 0;

			{ node n;
				forall_nodes(n,the_graph) {
					old_new_map[number[n]] = new_number[n];
//					std::cout << "old2new " << number[n] << " " << new_number[n] << std::endl;
				}
			}
    	}

        unsigned int* new_lp1 = new unsigned int[number_of_nodes+1];
        unsigned int* new_lp2 = new unsigned int[number_of_nodes+1];
        unsigned int* new_father = new unsigned int[number_of_nodes+1];

        { node n;
			forall_nodes(n,the_graph) {
				const unsigned int v = number[n];
				const unsigned int new_v = old_new_map[v];
				new_lp1[new_v] = old_new_map[lowpoint_one[v]];
				new_lp2[new_v] = old_new_map[lowpoint_two[v]];
				new_father[new_v] = old_new_map[father[v]];
			}
        }

        delete[] lowpoint_one;
        delete[] lowpoint_two;
        delete[] father;
        lowpoint_one = new_lp1;
        lowpoint_two = new_lp2;
        father = new_father;

        { node n;
			forall_nodes(n,the_graph) {
				number[n] = old_new_map[number[n]];
			}
        }

        delete[] old_new_map;


#ifdef DOTS
        std::fstream f("./step_2.dot", std::ios::trunc | std::ios::out);
        to_dot(f);
        f.close();
#endif
    }





    /* Step 1 in Hopcroft Tarjan algorithm. Compute the graph information.
     * Things changed:
     * - flag switched
     * - instead of looking at adjacent nodes I traverse adjacent edges, to be able to mark edges as fronds
     *
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
    unsigned int dfs_one(const node& node_v, const unsigned int parent, unsigned int next_number, bool* const seen_edge_to_parent) {
#ifdef DFSCOUT
		std::cout << "DFS " << node_v->id() << " " << node_v << " ";
#endif
		assert(next_number>0);
        number[node_v] = next_number++;

        assert(number[node_v]>0);
        assert(number[node_v]<=number_of_nodes);
        const unsigned int v = number[node_v];
        node_at[v] = node_v;

        assert(parent < v);
        assert(v>0);

        // We initialize lowpoint one and two with v, in case no lower node can be reached
        lowpoint_one[v] = v;
        lowpoint_two[v] = v;

        // We know that v has itself in its subtree
        number_descendants[node_v] = 1;

        {   edge vw;
            forall_adj_edges(vw,node_v) {

                node node_w = opposite(vw,node_v);
#ifdef DFSCOUT
                std::cout << "DFS edge "<< node_v->id() << " -> " << node_w->id() << std::endl;
#endif


                unsigned int w = number[node_w];

                if (w==0) {
                    //if w==0 we haven't seen w before, so we do a dfs on w
                    next_number = dfs_one(node_w, v,next_number,seen_edge_to_parent);

                    //during the dfs w gets assigned a number
                    w = number[node_w];
                    assert(w>0);

                    father[w] = v;
#ifdef DFSCOUT
                    std:: cout << "father " << w-1 << " = " << v-1 << std::endl;
#endif

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

                    /* Here I check for biconnectedness. If a graph is biconnected, v -> w => lowpoint_one[w] < v, unless v=1. In that case lowpoint_one[w] = v = 1
                     * If lowpoint_one[w] >= v, we can cut v out of the graph and separate w (and its descendants) from the rest, since there is no path back beyond v
                     * This biconnectedness test is most likely not exhaustive. I'm waiting for a bugfix in LEDA TODO
                     */
                    if (is_biconnected) {
    					if (v!=1) {
    						is_biconnected &= lowpoint_one[w] < v;
    					} else {
    						is_biconnected &= lowpoint_one[w] == v;
    					}
    					if (!is_biconnected) {
    						std::cout << "graph is not biconnected ";
    						if (v==1) {
    							std::cout << "v is root and" << w << " has lowpoint " << lowpoint_one[w] << std::endl;
    						} else {
    							std::cout << "edge from " << v  << " to " << w << " lowpoint_one[w] " << lowpoint_one[w] << std::endl;
    						}
    						articulation_point = node_at[v];
    					}
                    }

                } else if (w < v && ((w!=parent) || seen_edge_to_parent[v]) ) {
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

                    is_frond[vw] = true;
#ifdef DFSCOUT
                    std :: cout << "\t frond! " << (w<v) << " " << (w!=parent) << " " << seen_edge_to_parent[v] << std::endl;
#endif

                    /*
                     * A frond allows us to go back up in the graph, so we need to update the lowpoints
                     */
                    if (w<lowpoint_one[v]) {
                        lowpoint_two[v] = lowpoint_one[v];
                        lowpoint_one[v] = w;
                    } else if (w>lowpoint_one[v]) {
                        lowpoint_two[v] = min(lowpoint_two[v],w);
                    }

                } else {
#ifdef DFSCOUT
                	std::cout << "\t already seen " << std::endl;
#endif
                }


                // If we see another edge to u, it's a frond, as there is a path .
                if (w==parent) {
#ifdef DFSCOUT
                	std::cout << "DFS: seen_edge_to_parent "  << node_at[parent]->id() <<  std::endl;
#endif
                	seen_edge_to_parent[v] = true;
                }


            }
        }

        return next_number;
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
    void pathfinder(node node_v, bool& at_path_start, unsigned int nodes_to_be_seen, node_array<unsigned int>& new_number) {
        assert(nodes_to_be_seen <=number_of_nodes);

        /* because we want the numbers in reverse, but need them before we come back (for the highpoints) we count downwards.
         * so we decrease the number each time we discover a vertex. Since we have number_descendants[v]-1 vertices to
         * see when we first visit v, we will decrease nodes_to_be_seen by that amount.
         */
        assert(nodes_to_be_seen>=(number_descendants[node_v]-1));
        const unsigned int v = nodes_to_be_seen - (number_descendants[node_v] - 1);
        assert(v<=number_of_nodes);
        new_number[node_v] = v;
        assert(new_number[node_v]>0);

        node_at[v] = node_v;
#ifdef PATHFINDER_COUT
        std::cout << "Pathfinder " << node_v->id() << " Newnum " << v << " to be seen " << nodes_to_be_seen << " nd " << number_descendants[node_v] << std::endl;
#endif


        {   edge vw;
            forall_adj_edges(vw,node_v) {
                /* In this loop we check whether the edge vw is the first on one of the paths */

                const node node_w = opposite(vw,node_v);
                const unsigned int w = new_number[node_w];

                if ((w==0 || w>v)&& is_frond[vw]) {
#ifdef PATHFINDER_COUT
                	std::cout << "Skip frond in wrong direction " << std::endl;
#endif
                	continue;
                }

                if (w>0 && w<v && !is_frond[vw]) {
#ifdef PATHFINDER_COUT
                	std::cout << "Skip edge to parent " << std::endl;
#endif
                	continue;
                }

				if (at_path_start) {
#ifdef PATHFINDER_PATH_COUT
					std::cout << std::endl << "Path ";
#endif
					is_first_on_path[vw] = at_path_start; //if we see edges again that are true, don't set it to false
					at_path_start=false;
				}
#ifdef PATHFINDER_PATH_COUT
				std::cout << node_v->id() << "->" << node_w->id() << " " ;
#endif

                if (!is_frond[vw] ) {
#ifdef PATHFINDER_COUT
                		std::cout << "\ttree edge " << std::endl;
#endif
                		pathfinder(node_w, at_path_start, nodes_to_be_seen, new_number);
                		assert(new_number[node_w]>0);
						nodes_to_be_seen -= number_descendants[node_w];
                } else {
                    /* if we encounter a frond, it is the last edge on the path, so set s to true to mark the next edge as a start of a path.
                     * we also update the highpoint. We don't need to check whether it is smaller than before because of the ordering in which we traverse the graph
                     */

#ifdef PATHFINDER_COUT
                	std::cout << "\tfrond " << node_v->id() << " " << node_w->id() <<  std::endl;
#endif
                    if (highpoint[node_w] == 0) {
#ifdef PATHFINDER_COUT
                    	std::cout << "Frond from " << node_v->id() << " (" << v << ")" << " to " << node_w->id() << " (" << w << "): update highpoint " << w << " to " << v << std::endl;
#endif
                        highpoint[node_w] = v;
                    }
#ifdef PATHFINDER_COUT
                    else {
                    	std::cout << "Frond from " <<  node_v->id() << " to " << node_w->id() << ": highpoint stays at " << highpoint[node_w] <<  std::endl;
                    }
#endif

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
#ifdef PATHSEARCH_COUT
        std::cout<< "Pathsearch: " << v << " " << node_at[v]->id() << std::endl;
#endif

        //unsigned int num_unexplored_adj_edges = the_graph.degree(node_v);

        edge vw;
        forall_adj_edges(vw,node_v) {
            //num_unexplored_adj_edges--;
            const node node_w = opposite(vw,node_v);
            const bool frond = is_frond[vw];
            const unsigned int w = number[node_w];
#ifdef PATHSEARCH_COUT
            std::cout << "\tEdge " << source(vw)->id() << " " << target(vw)->id() << std::endl;
#endif

            if (w < v && !frond) {
#ifdef PATHSEARCH_COUT
				std::cout << "Skip backedge to parent " << source(vw)->id() << " " << target(vw) -> id() << " from " << node_v->id() << std::endl;
#endif
				continue;
            }

            if (is_frond[vw] && v < w) {
#ifdef PATHSEARCH_COUT
            	std::cout << "Skip frond in the wrong direction " << node_v->id() << " " << node_w->id() << std::endl;
#endif
				continue;
            }


            if (is_first_on_path[vw])
                update_stack(vw,node_w, w, v,current_t_stack);

            if (frond) {
#ifdef PATHSEARCH_COUT
            	std::cout << "Skip frond " << source(vw)->id() << " " << target(vw) -> id() << std::endl;
#endif
            	continue;
            }

            if (!pathsearch(node_w,s1,s2))
                return false;
#ifdef PATHSEARCH_COUT
            std::cout << "Continuing with edge " << source(vw)->id() << " " << target(vw)->id() << std::endl;
#endif

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
                {   const node node_child_w = the_graph.adj_nodes(node_w).contents(the_graph.adj_nodes(node_w).first());
                    // it should also be ok to check whether the first out-edge of w is a frond.
                    // I don't know why we need to check number[node_child_w] > w. TODO
                    const bool deg_two_case =  (the_graph.degree(node_w) == 2 && number[node_child_w] > w);

                    if (deg_two_case) {
#ifdef PATHSEARCH_COUT_SEPPAIR
                    	std::cout << v << " and " << number[node_child_w] << " are a type-2 pair, " << w << " has degree 2" <<std::endl;
#endif
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
#ifdef PATHSEARCH_COUT
                    	std::cout << "Not a type 2 pair (f[b]=v) (" <<  current_t_stack->top().first() << ", " << current_t_stack->top().second() << ", " << current_t_stack->top().third() << ")" << std::endl;
#endif

                        current_t_stack->pop();
                    } else {
#ifdef PATHSEARCH_COUT_SEPPAIR
                    	std::cout << number[node_v] << " and " << b << " are a type-2 pair" << std:: endl;
#endif
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
             */

            // i satisfied by choice of edge v -> w
            // if lowpoint_two[v] >= v and lowpoint_one[w] >= v removal of node v should disconnect the graph => not biconnected.
            {
            	const unsigned int b = v;
            	const unsigned int a = lowpoint_one[w];
#ifdef PATHSEARCH_COUT
            	std::cout << "Possible type 1 pair: " << a << " " << b << std::endl;
#endif

				if (/* iii */ lowpoint_two[w] >= b && a < b && (father[b] != 1 || number_descendants[node_v] - number_descendants[node_w] > 1))
				{
					s1 = node_at[a]; // ii
					s2 = node_v;
#ifdef PATHSEARCH_COUT_SEPPAIR
					std::cout << a << " and " << b << " are a type-1 pair (";
					std::cout << s1->id() << ", " << s2->id() << ")" << std::endl;
					std::cout << (lowpoint_two[w] >= a) << " " << (a < b) << " " << (father[b]!=1) << " " << (number_descendants[node_v] - number_descendants[node_w]) << std::endl;
					std::cout << "w = " << w << std::endl;
#endif
					return false;
				}
            }

            if (is_first_on_path[vw]) {
#ifdef PATHSEARCH_COUT
            	std::cout << "Deleting stack ";
            	print_stack(*current_t_stack);
#endif
            	delete current_t_stack;
                current_t_stack = t_stacks.pop();
            }

            /*
             * This is another step for maintaining the t_stack. When we backtrack over the tree edge v -- w it is possible that (v)
             *
             * v) and every frond x-->v with a<v<b and b -> w -*>x has lowpoint_one[w] >= a
             *
             * becomes violated. Remember
             *
             * highpoint[v], for nodes at which fronds end, the node with the lowest number at the source of a frond.
             *
             * So if highpoint[v] > h for a potential sep. pair (h,a,b), a<v<b, that means that in the corresponding split
             * component there is a frond back into the part between a and b and hence it can't be a valid separation pair.
             */
#ifdef PATHSEARCH_COUT
            std::cout << "Deleting things from stack with h < highpoint["<< v << "] = " << highpoint[node_v] << std::endl;
#endif
            while (!current_t_stack->empty()
                    && current_t_stack->top().second() != v
                    && current_t_stack->top().third() != v
                    && highpoint[node_v] > current_t_stack->top().first()) {
#ifdef PATHSEARCH_COUT
				std::cout << "Not a type 2 pair (hp) (" <<  current_t_stack->top().first() << ", " << current_t_stack->top().second() << ", " << current_t_stack->top().third() << ")" << std::endl;
#endif
                assert(current_t_stack->top().second() < v);
                assert(current_t_stack->top().third() > v);
                current_t_stack->pop();
            }


        }

        return true;
    }

    void print_stack(t_stack_type stack) {
    	while(!stack.empty()) {
    		std::cout << "(" << stack.top().first() << ", " << stack.top().second() << ", " << stack.top().third() << ")";
    		stack.pop();
    	}
    	std::cout << std::endl;
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
    void update_stack(const edge e, const node node_w, const unsigned int w, const unsigned int v, t_stack_type*& current_t_stack) {
#ifdef STACK_COUT
    	std::cout << "Stackupdate " << current_t_stack->size() << " " << t_stacks.size() << " edge " << v << " " << w << std::endl;
    	std::cout << "Stack ";
    	print_stack(*current_t_stack);
#endif
        if (!is_frond[e]) {
        	assert(v<w);

            /* When we start a new path with the tree edge v -> w it's possible to violate (iv) for possible pairs (a,b) on the stack. So we check that lowpoint_one[w] <= a and pop,
             * if necessary.
             *
             * If nothing was popped (lowpoint_one[w],v) is a possible separation pair (why? TODO). We push it on the stack, together with the highest node number in w's subtree.
             * Otherwise, if some triple (h,a,b) was deleted last, we get a new possible sep. pair with (lowpoint_one[w],b), as this pair obviously satisfies (iv). We push that
             */

            // in case nothing gets deleted, we want to push these values:

            unsigned int highest_vertex_number = w + number_descendants[node_w] -1;
            three_tuple<unsigned int, unsigned int, unsigned int> hab(highest_vertex_number, lowpoint_one[w], v);

            // remove all tuples from the stack that violate (iv), remember the last one
            three_tuple<unsigned int, unsigned int, unsigned int> last_popped = hab;
            while(!current_t_stack->empty() && current_t_stack->top().second() > lowpoint_one[w]) {
                last_popped = current_t_stack->pop();
                highest_vertex_number = max(highest_vertex_number, last_popped.first());
#ifdef STACK_COUT
                std::cout << "\tPop " << last_popped.first() << " " << last_popped.second() << " " << last_popped.third() << std::endl;
#endif
            }

            // (lowpoint_one[w],b) is a possible sep. pair
            hab.first() = max(highest_vertex_number,last_popped.first());
            hab.third() = last_popped.third();

#ifdef STACK_COUT
            std::cout << "Possible type 2 pair: " << hab.first() << " " << hab.second() << " " << hab.third() << std::endl;
#endif
            assert(current_t_stack->empty() || (current_t_stack->top().second() <= hab.second() && current_t_stack->top().third() >= hab.third()));
            assert(hab.second() <= v);
            assert(hab.third() >= v);
            assert(hab.second()<=hab.third());
            current_t_stack->push(hab);

            // Since we will continue with examining things on the new cycle, we save this t_stack for later, when we return.
            t_stacks.push(current_t_stack);

            current_t_stack = new t_stack_type();

        } else {
        	assert(v>w);
            /* In the case of a frond, the new path consists of just one edge, so the top pair of the stack is ok, if w >= a (v). So we pop all triples with w<a.
             *
             * If nothing was deleted (w,v) is a possible type-2 pair. Since v is processed last in this split component it has the highest number and we push (v,w,v)
             * Otherwise let (h,a,b) be the last triple that got deleted. Then (w,b) is a possible new pair. The highest number for that pair is the maximum that was deleted.
             */

            unsigned int highest_vertex_number = v;
            three_tuple<unsigned int, unsigned int, unsigned int> hab(v,w,v);

            three_tuple<unsigned int, unsigned int, unsigned int> last_popped = hab;
            //pop all triples (h,a,b) with a>w, remember the last one
            while(!current_t_stack->empty() && w < current_t_stack->top().second()) {
                last_popped = current_t_stack->pop();
                highest_vertex_number = max(highest_vertex_number, last_popped.first());
#ifdef STACK_COUT
                std::cout << "\tPop " << last_popped.first() << " " << last_popped.second() << " " << last_popped.third() << std::endl;
#endif
            }

            hab.first() = highest_vertex_number;
            hab.third() = last_popped.third();
#ifdef STACK_COUT
            std::cout << "Possible type 2 pair: " << hab.first() << " " << hab.second() << " " << hab.third() << std::endl;
#endif
            assert(current_t_stack->empty() || (current_t_stack->top().second() <= hab.second() && current_t_stack->top().third() >= hab.third()));
            assert(hab.second() <= v);
            assert(hab.third() >= v);
            assert(hab.second()<=hab.third());
            current_t_stack->push(hab);
            assert(hab.first()>=v);

            // Since this segment is done processing after this edge, we don't have to start a new t_stack

        }
    }
};

bool hopcroft_tarjan_is_triconnected(const ugraph& g, node& s1, node& s2) {
    if (!is_connected(g))
        return false;

    ugraph H;
    CopyGraph(H,g);
    palm_tree p(H);
    return p.is_triconnected(s1,s2);
}

bool hopcroft_tarjan_is_triconnected(const ugraph& g) {
    node s1,s2;
    return hopcroft_tarjan_is_triconnected(g,s1,s2);
}

bool hopcroft_tarjan_is_triconnected_nc(ugraph& g, node& s1, node& s2) {
    if (!is_connected(g))
        return false;

    palm_tree p(g);
    return p.is_triconnected(s1,s2);
}



void test(ugraph& g) {
    palm_tree p(g);
    std::fstream f("./the_palm_tree.dot", std::ios::trunc | std::ios::out);
    p.to_dot(f);
}
