#include "triconnectivity.hpp"
#include "dfs.hpp"

using namespace leda;

/* A graph is triconnected if there is no pair of nodes u,v s.t.
 * the graph is no longer connected if u and v are removed. i.e.
 * there exists a pair of nodes x,y s.t. all paths between x,y pass through u or v.
 */
bool naive_is_triconnected(ugraph& g) {

	if (!is_connected(g))
		return false;
	node u;
	forall_nodes(u,g) {
		int degree = g.degree(u);
		list<edge> adj_edges = g.adj_edges(u);
		g.hide_node(u);

		node v;
		forall_nodes(v,g) {
			int degree = g.degree(v);
			list<edge> adj_edges = g.adj_edges(v);
			g.hide_node(v);

			if (!is_connected(g)) { //u,v is a separation pair
				g.restore_all_nodes();
				g.restore_all_edges();
				return false;
			}

			g.restore_node(v);
			g.restore_edges(adj_edges);
			assert(!g.is_hidden(v));
			assert(degree==g.degree(v));

		}
		g.restore_node(u);
		g.restore_edges(adj_edges);
		assert(degree==g.degree(u));

	}
	return true;
}
