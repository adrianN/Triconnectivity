#include "certificate.hpp"
#include "chain_edge_iterator.hpp"
#include "utilities.hpp"


certificate::certificate(ugraph  & graph,   schmidt_triconnectivity* d) :
	the_graph(graph),
	decomposition(d),
	created_by_chain(new_graph,-1),
	orig_2_new(the_graph,NULL),
	new_2_orig(new_graph,NULL)
	{}

certificate::~certificate() {
}

bool certificate::add_bg_path(list<edge> const & edges) {
	std::cout << "\t****** BG PATH: " ;
	node first_node = NULL;
	node last_node = NULL;
	node current_node = NULL;
	edge e;

	chains.push_back(new list<edge>());

	unsigned int prev_created_nodes = 0;
	forall(e,edges) {
		std::cout << decomposition->dfi(source(e)) << " -> " << decomposition->dfi(target(e)) << " ";
		if (current_node == NULL)
			current_node = source(e);
		else
			current_node = opposite(current_node,e);

		if (first_node == NULL)
			first_node = current_node;

		if (last_node == NULL) {
			last_node = current_node;
			continue;
		} else {
			node new_node;
			if (orig_2_new[current_node] == NULL) {
				orig_2_new[current_node] = new_node = new_graph.new_node();
				new_2_orig[new_node] = current_node;
			} else {
				prev_created_nodes++;
			}
			node prev_node = orig_2_new[last_node];
			assert(prev_node !=NULL);
			edge added_edge = new_graph.new_edge(prev_node,new_node);
			list<edge>::item it = chains.back()->append(added_edge);
			le_edges[added_edge] = new pair<list<edge>*,list<edge>::item>(chains.back(), it);
		}
	}

	assert(prev_created_nodes == 2);
	assert(first_node != last_node);

	endvertices.push_back(new pair<node,node>(first_node,last_node));

	std::cout << std::endl;


	return true;
}

bool certificate::add_bg_path(const chain * c) {
	list<edge> l;
	chain_edge_iterator it(c, decomposition);
	for(edge e = it.next(); e!=NULL; e=it.next()) {
		l.append(e);
	}
	add_bg_path(l);
	return true;
}

bool certificate::verify() {
	if (!graphs_isomorphic(the_graph, new_graph, new_2_orig)) {
		return false;
	}
	for(int i = chains.size()-1; i>=0; i--) {
		assert(chains[i]->size() == 1); // there is one edge left in this chain.
		const edge to_be_deleted = chains[i]->head();

		const node a = endvertices[i]->first;
		const node b = endvertices[i]->second;
		delete endvertices[i];
		endvertices[i] = NULL;

		const bool a_smooth = new_graph.degree(a) == 2;
		const bool b_smooth = new_graph.degree(b) == 2;
		unsigned int type = 0;
		if (a_smooth) type++;
		if (b_smooth) type++;

		//check BG conditions

		switch(type) {
		case 2: {		//subdivide 2, connect

			node a_neighbours[2];
			node b_neighbours[2];
			{	edge e;
				unsigned int i = 0;
				forall_adj_edges(e,a) {
					a_neighbours[i++] = opposite(e,a);
				}
				i=0;
				forall_adj_edges(e,b) {
					b_neighbours[i++] = opposite(e,b);
				}
			}
			if (a_neighbours[0] < a_neighbours[1]) {
				std::swap(a_neighbours[0], a_neighbours[1]);
			}
			if (b_neighbours[0] < b_neighbours[1]) {
				std::swap(b_neighbours[0], b_neighbours[1]);
			}
			if (a_neighbours[0] == b_neighbours[0] && a_neighbours[1] == b_neighbours[1]) {
				assert(false && "chain between parallel links");
			}
			// update le_edges and the chains
			//first remove the smoothened edges from the chains
			{ 	list<edge>* lists[2][2];
				unsigned int i=0, j=0;

				{
					edge e;
					forall_adj_edges(e,a) {
						pair<list<edge>*, list<edge>::item>* p = le_edges[e];
						p->first->del_item(p->second);
						lists[j][i++] = p->first;
					}
					j++;
					forall_adj_edges(e,b) {
						pair<list<edge>*, list<edge>::item>* p = le_edges[e];
						p->first->del_item(p->second);
						lists[j][i++] = p->first;
					}
				}
				// then remove them from the graph
				edge es[2];
				es[0] = smoothe(new_graph,a);
				es[1] = smoothe(new_graph,b);

				//and insert the new edges into the chains
				for(unsigned int i = 0; i<2; i++) {
					assert(lists[i][1] == lists[i][0]);
					list<edge>::item it = lists[i][0]->append(es[i]);
					le_edges[es[i]] = new pair<list<edge>*,list<edge>::item>(lists[i][0],it);
				}
			}

		} break;
		case 1: {		//subdivide 1, connect
			node third, sub;
			if (a_smooth) {
				third = b;
				sub = a;
			} else {
				assert(b_smooth);
				third = a;
				sub = b;
			}
			node neighbours[2];

			{	edge e;
				unsigned int i = 0;
				forall_adj_edges(e,sub) {
					neighbours[i++] = opposite(e,sub);
				}
			}
			if (neighbours[0] == third || neighbours[1] == third) {
				assert(false && "third node is one of the endpoints");
			}

			//update le_edges and the chains
			//first delete smoothened edges from the chain
			{ 	list<edge>* lists[2];
				unsigned int i = 0;
				{	edge e;
					forall_adj_edges(e,sub) {
						pair<list<edge>*, list<edge>::item>* p = le_edges[e];
						p->first->del_item(p->second);
						lists[i++] = p->first;
					}
				}
				//then from the graph
				edge e1 = smoothe(new_graph,sub);

				//and insert the new edge into the chain
				assert(lists[0] == lists[1]);
				list<edge>::item it = lists[0]->append(e1);
				le_edges[e1] = new pair<list<edge>*, list<edge>::item>(lists[0],it);
			}
		} break;
		case 0: {		//connect two is always okay, because of the checks we performed before adding this chain.
		} break;
		default:
			assert(false);
		}

		// delete the chain
		delete(le_edges[to_be_deleted]);
		le_edges[to_be_deleted] = NULL;

		new_graph.del_edge(to_be_deleted);

		delete(chains[i]);
		chains[i] = NULL;

	}

	return true;
}
