#include "triconnectivity.hpp"
#include "dfs.hpp"

using std::auto_ptr;

using namespace leda;

/* A graph is triconnected if there is no pair of nodes u,v s.t.
 * the graph is no longer connected if u and v are removed. i.e.
 * there exists a pair of nodes x,y s.t. all paths between x,y pass through u or v.
 */
bool naive_is_triconnected(ugraph& g) {
    if (!is_connected(g)) return false;
    auto_ptr<list<two_tuple<node,node> > > sep_pairs =naive_separation_pairs(g);
    bool b = sep_pairs->length()==0;

    return b;
}


/* Return a list of all separation pairs of g. Precondition: g is connected */
auto_ptr<list<two_tuple<node,node> > > naive_separation_pairs(ugraph& g) {

    assert(is_connected(g));
    node u;
    auto_ptr<list<two_tuple<node,node> > > pairs(new list<two_tuple<node,node> >());

    node* nodes_of_h = new node[g.number_of_nodes()];
    { int i = 0;
    forall_nodes(u,g) {
    	nodes_of_h[i++] = u;
    }
    }

    const int FUCKING_NUMBER_OF_FUCKING_NODES = g.number_of_nodes();
    for(int i=0; i<FUCKING_NUMBER_OF_FUCKING_NODES; i++) {
    	node u = nodes_of_h[i];
    	list<edge> adj_edges = g.adj_edges(u);

    	g.hide_node(u);

    	for(int j=i+1; j<FUCKING_NUMBER_OF_FUCKING_NODES; j++) {
    		node v = nodes_of_h[j];
    		list<edge> adj_edges = g.adj_edges(v);

			g.hide_node(v);
//			std::cout << "naive: " << u->id() << " " << v->id() << " " << i << " " << j << "\n";

            if (!is_connected(g)) { //u,v is a separation pair
                pairs->append(two_tuple<node,node>(u,v));
            }

            g.restore_node(v);
            g.restore_edges(adj_edges);
    	}

    	g.restore_node(u);
        g.restore_edges(adj_edges);
    }

    delete[] nodes_of_h;
    return pairs;
}
