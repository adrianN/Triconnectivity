#include "triconnectivity.hpp"
#include "LEDA/core/slist.h"
#include <ostream>
#include <vector>
#include <memory>

using namespace leda;
using std::endl;
using std::istream;
using std::ostream;
using std::auto_ptr;
using std::cout;
using std::endl;

//#define DETERMINISTIC_GLUE

void make_simple(ugraph& u) {
	node n;
	forall_nodes(n,u) {
		node_array<bool> seen(u,false);
		edge e;
		forall_adj_edges(e,n) {
			if (seen[opposite(e,n)])
				u.hide_edge(e);
			else
				seen[opposite(e,n)] = true;
		}
	}
}

edge smoothen(ugraph & g, node n) {
	edge e=NULL;
	assert(g.degree(n)==2);
	node neighbours[2]={NULL,NULL};
	{	unsigned int i = 0;
		edge e;
		forall_adj_edges(e,n) {
			neighbours[i++] = opposite(e,n);
		}
	}

	e = g.new_edge(neighbours[0], neighbours[1]);

	g.del_node(n);

	return e;
}

//this may or may not depend on undefined behaviour.
int edge_to_int(const node n1, const node n2) {
	int n = 0;
	assert(n1->id() < USHRT_MAX);
	assert(n2->id() < USHRT_MAX);
	unsigned short id1 = n1->id();
	unsigned short id2 = n2->id();
	n = id1;
	n *= (USHRT_MAX+1);
	n += id2;
	return n;
}

void simple_to_dot(const ugraph&,ostream&);

bool graphs_isomorphic(ugraph  & g1, ugraph  & g2, node_array<node> const & map_2_to_1) {
	if (g1.number_of_edges() != g2.number_of_edges() || g1.number_of_nodes() != g2.number_of_nodes())
		return false;

	edge_array<int> g1_order(g1);
	edge_array<int> g2_order(g2);
	{ 	edge e;
		forall_edges(e,g1) {
			if (source(e)->id() > target(e)->id())
				g1.rev_edge(e);
			g1_order[e] = edge_to_int(source(e),target(e));
		}
		forall_edges(e,g2) {
			if (map_2_to_1[source(e)] > map_2_to_1[target(e)])
				g2.rev_edge(e);
			g2_order[e] = edge_to_int(map_2_to_1[source(e)], map_2_to_1[target(e)]);
		}
	}
	//g1.bucket_sort_edges(g1_order);
	//g2.bucket_sort_edges(g2_order);
	g1.sort_edges(g1_order);
	g2.sort_edges(g2_order);
	edge e1,e2;
	e1 = g1.first_edge();
	e2 = g2.first_edge();
	while((e1 !=NULL && e2 != NULL) && (source(e1) == map_2_to_1[source(e2)] && target(e1) == map_2_to_1[target(e2)])) {

		e1 = g1.succ_edge(e1);
		e2 = g2.succ_edge(e2);
	}
#ifndef NDEBUG
	if (e1!=NULL || e2!=NULL) {
		std::fstream f1("./g1.dot", std::ios::out);
		std::fstream f2("./g2.dot",std::ios::out);
		simple_to_dot(g1,f1);
		simple_to_dot(g2,f2);
		if (e1!=NULL && e2 !=NULL)
		std::cout << "Edge " << source(e1) ->id() << ", " << target(e1)->id() << " and " << source(e2)->id() << "," << target(e2)->id() << "(" << map_2_to_1[source(e2)]->id() << "," << map_2_to_1[target(e2)]->id() << ")" << std::endl;
	}
#endif
	return (e1 == NULL && e2 == NULL);
}

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

	unsigned int num_nodes;
	assert(input.good());
	assert(!input.eof());
	num_nodes = input.get();
	if (num_nodes == '>' && input.peek() == '>') { //file starts with >>planar_code<<
		input.ignore(14);
		num_nodes = input.get();
		assert(num_nodes != '<');
	}

#ifdef GRAPHREAD_COUT
	std::cout << "\t number of nodes: " << num_nodes << std::endl;
#endif

	node* nodes = new node[num_nodes];

	if (num_nodes != (unsigned int)graph.number_of_nodes()) {
		graph.clear();
		for(unsigned int i=0; i<num_nodes; i++) {
			nodes[i] = graph.new_node();
		}
	} else {
		node n;
		unsigned int i=0;
		forall_nodes(n,graph) {
			nodes[i++] = n;
		}
		graph.del_all_edges();
	}


	for(unsigned int i=0; i<num_nodes; i++) {
		unsigned int adj_node;
		adj_node = input.get();
		while (adj_node!=0) {
			if (adj_node-1>i) {
#ifdef GRAPHREAD_COUT
				std::cout<< "\tEdge from " << i << " to " << (int)adj_node-1 << std::endl;
#endif

				graph.new_edge(nodes[i],nodes[adj_node-1]);

			}
			adj_node = input.get();
			assert(input.good());
		}
	}
	delete[] nodes;
	input.peek();
	return input;
}

node merge_nodes(ugraph& g, slist<node> nodes) {
	node new_node = g.new_node();
	slist<node> edge_to;
	node_array<bool> no_more_edges_to(g,false);

	//collect all edges, ban nodes in the list
	{	node u;
		forall(u, nodes) {
			no_more_edges_to[u] = true;
			node v;
			forall_adj_nodes(v,u) {
				edge_to.append(v);
			}
		}
	}

	{//iterate over edges, add new edges to new_node
		node v;
		forall(v,edge_to) {
			if (no_more_edges_to[v]) continue;

			no_more_edges_to[v] = true;

			g.new_edge(new_node,v);
		}
	}

	{	node v;
		forall(v,nodes) {
			//g.del_node(v);
			g.hide_node(v);
		}
	}
	return new_node;
}

node merge_nodes(ugraph& g, node one, node two) {
	node_array<bool> edge_to(g,false);
	node min, max;
	if (g.degree(one) > g.degree(two)) {
		min = two;
		max = one;
	} else {
		min = one;
		max = two;
	}
	edge e;
	forall_adj_edges(e,max) {
		node opp = opposite(e,max);
		edge_to[opp] = true;
	}

	forall_adj_edges(e,min) {
		node opp = opposite(e,min);
		if (!edge_to[opp] && max != opp) {
			//g.new_edge(max,opp);
			g.move_edge(e,max,opp);
			edge_to[opp] =true;
		}
	}
	g.del_node(min);

	return max;
//	node new_node = g.new_node();
//	node_array<bool> edge_to(g,false);
//	edge e;
//	forall_adj_edges(e,one) {
//		node opp = opposite(e,one);
//		if (!edge_to[opp] && new_node!=opp) {
//			g.new_edge(new_node,opp);
//			edge_to[opp] = true;
//		}
//	}
//	forall_adj_edges(e,two) {
//		node opp = opposite(e,two);
//		if (!edge_to[opp] && new_node!=opp) {
//			g.new_edge(new_node,opp);
//			edge_to[opp] = true;
//		}
//	}
//	g.hide_node(one);
//	g.hide_node(two);
////	g.del_node(one);
////	g.del_node(two);
//
//	return new_node;
}

/* Takes two graphs and indentifies two nodes from each to form a new graph with exactly one separation pair. The result is stored in one. The separation pair is stored in the node arguments */
void glue_graphs(ugraph& one, ugraph& two, node& merge_one, node& merge_two) {
#ifndef NDEBUG
	const unsigned int number_one_nodes = one.number_of_nodes();
#endif
	const unsigned int number_two_nodes = two.number_of_nodes();

#ifdef DETERMINISTIC_GLUE
	const node point_one = one.first_node();
	node point_two = one.last_node();
	assert(point_one != point_two);
#else
	const node point_one = one.choose_node();
	node point_two;
	do {
		point_two = one.choose_node();
	} while(point_one == point_two);
#endif
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
#ifdef DETERMINISTIC_GLUE
	src_two.set_seed(42);
#endif

	unsigned int other_graph_point_one, other_graph_point_two;
	src_two >> other_graph_point_one;
	do {
		src_two >> other_graph_point_two;
	} while(other_graph_point_one == other_graph_point_two);

	merge_one = merge_nodes(one,point_one, new_nodes[other_graph_point_one]);
	merge_two = merge_nodes(one,point_two, new_nodes[other_graph_point_two]);

//	std::cout << "Merge nodes: " << merge_one->id() << " " << merge_two->id() << " from " << std::endl;
//	std::cout << point_one->id() << " " << new_nodes[other_graph_point_one]->id() << std::endl;
//	std::cout << point_two->id() << " " << new_nodes[other_graph_point_two]->id() << std::endl;

	delete[] new_nodes;

	assert((unsigned int)one.number_of_nodes() == number_one_nodes + number_two_nodes - 2);
}

void write_planar_code(ugraph& g, std::ostream& out) {
	node_array<unsigned char> ids(g,0);
	unsigned char id=1;
	node n;
	forall_nodes(n,g) {
		ids[n] = id++;
	}
	out << (unsigned char)g.number_of_nodes();
	forall_nodes(n,g) {
//		std::cout << ids[n] << std::endl;
		edge e;
		forall_adj_edges(e,n) {
//			std::cout << "\t" << ids[opposite(e,n)] << std::endl;
			out << ids[opposite(e,n)];
		}
		out << '\0';
	}
}

auto_ptr<ugraph> triconnected_graph(const unsigned int n) {
	auto_ptr<ugraph> g(new ugraph());
	assert(n>=4);

	//start with k4
	node* nodes = new node[n];
	unsigned int node_ptr=4;

	{
		const unsigned int num_nodes = 4;
		const unsigned int max_degree = 3;
	    for(unsigned int i=0; i<num_nodes; i++) {
	        nodes[i] = g->new_node();
	    }
	    int edges[num_nodes][max_degree] = {
	    		{1,3,2},
	    		{2,3,-1},
	    		{3,-1,-1},
	    		{-1,-1,-1}
	    };

	    for(unsigned int i=0; i<num_nodes; i++) {
	        for(unsigned int j=0; j<max_degree; j++) {
	            if (edges[i][j]>=0) {
	                g->new_edge(nodes[i],nodes[edges[i][j]]);
	            }
	        }
	    }
	    {   edge_array<int> edges(*g);
	        edge e;
	        forall_edges(e,*g) {
	            edges[e] = target(e)->id();
	        }
	        g->bucket_sort_edges(edges);
	    }
	}

	node neighbours[3];
	random_source r;
	for(unsigned int i = 0; i<n-4; i++) {
		if (i%(n/20)==0) {
			std::cout << '-';
			std::cout.flush();
		}
		unsigned int u,v,w;

		do {
			u = r.get() % node_ptr;
			v = r.get() % node_ptr;
			w = r.get() % node_ptr;

//			neighbours[0] = g->choose_node();
//			neighbours[1] = g->choose_node();
//			neighbours[2] = g->choose_node();
//			assert(neighbours[0]!=0);
//			assert(neighbours[1]!=0);
//			assert(neighbours[2]!=0);
		} while (u==v || v==w || w == u);

		neighbours[0] = nodes[u];
		neighbours[1] = nodes[v];
		neighbours[2] = nodes[w];

		nodes[node_ptr] = g->new_node();
		g->new_edge(nodes[node_ptr],neighbours[0]);
		g->new_edge(nodes[node_ptr],neighbours[1]);
		g->new_edge(nodes[node_ptr],neighbours[2]);
		node_ptr++;
	}
	delete[] nodes;
	return g;
}

void reduce(ugraph& g) {
	char c=0;
	ugraph h(g);
	char i='a';
	while (c!='x') {
		node s1=NULL,s2=NULL;
		bool tc = hopcroft_tarjan_is_triconnected_nc(g,s1,s2);
		std::cout << "tc " << tc << std::endl;
		if (!tc) {
			std::cout << s1->id() << std::endl;
			std::cout << s2->id() << std::endl;
		}

		i++;
		{
			std::string s = "reduced";
			s+=i;
			s+=".tri";
			std::fstream reduced(s.c_str(),std::ios::out);
			write_planar_code(g,reduced);
			reduced.close();
		}
		std::cin >> c;
		switch(c) {
		case 'd': {
			int n;
			std::cin >> n;
			node no;
			CopyGraph(h,g);
			forall_nodes(no,g) {
				if (no->id() == n) {
					g.del_node(no);
					break;
				}
			}
			break;
		}
		case 'm': {
			int n,m;
			std::cin >> n;
			std::cin >> m;
			node non=0,nom=0;
			node no;
			forall_nodes(no,g) {
				if (no->id() == n)
					non = no;
				if (no->id() == m)
					nom = no;
			}
			CopyGraph(h,g);

			merge_nodes(g,non,nom);
			break;
		}
		case 'e': {
			int n,m;
			std::cin >> n;
			std::cin >> m;
			edge e;
			edge f=0;
			forall_edges(e,g) {
				if ((source(e)->id() == n && target(e)->id() == m) || (source(e)->id() == m && target(e)->id() == n)) {
					f = e;
					break;
				}
			}
			CopyGraph(h,g);

			g.del_edge(f);
			break;
		}
		case 'u': {
			CopyGraph(g,h);
			break;
		}
		}

	}

}

