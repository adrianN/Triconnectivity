#include "certificate.hpp"
#include "chain_node_iterator.hpp"
#include "utilities.hpp"

//#define BGCOUT

certificate::certificate(ugraph  & graph,   schmidt_triconnectivity* d) :
	still_valid(true),
	the_graph(graph),
	decomposition(d),
	created_by_chain(new_graph,the_graph.number_of_nodes(), -1),
	le_edges(new_graph,the_graph.number_of_edges(), NULL),
	orig_2_new(the_graph,NULL),
	new_2_orig(new_graph,the_graph.number_of_nodes(),NULL)
	{}

certificate::~certificate() {
}

bool certificate::add_bg_path(list<node> const & nodes) throw() {
	if (!still_valid)
		return false;
#ifdef BGCOUT
	std::cout << "\t****** BG PATH: " ;
#endif
	node first_node = NULL;
	node last_node = NULL;
	node current_node = NULL;
	chains.push_back(new list<edge>());

	node n;
	unsigned int prev_created_nodes = 0;
	forall(n,nodes) {
#ifdef BGCOUT
		std::cout << decomposition->dfi(n) << " ";
#endif
		current_node = n;

		if (first_node == NULL)
			first_node = current_node;

		node new_node;
		if (orig_2_new[current_node] == NULL) {
			new_node = new_graph.new_node();
			orig_2_new[current_node] = new_node;
			new_2_orig[new_node] = current_node;
		} else {
			new_node = orig_2_new[current_node];
			prev_created_nodes++;
		}

		if (last_node !=NULL) {
			node prev_node = orig_2_new[last_node];
			assert(prev_node !=NULL);
			edge added_edge = new_graph.new_edge(prev_node,new_node);
			list<edge>::item it = chains.back()->append(added_edge);
			le_edges[added_edge] = new pair<list<edge>*,list<edge>::item>(chains.back(), it);
		}

		last_node = current_node;

	}

	// paths have only two nodes in common with the rest of the graph
	still_valid &= chains.size() == 1 || prev_created_nodes == 2;

	// bg paths mustn't form cycles
	still_valid &= first_node != last_node;

	endvertices.push_back(new pair<node,node>(first_node,last_node));
#ifdef BGCOUT
	std::cout << std::endl;
#endif

#ifndef NDEBUG
	if (!still_valid)
		std::cout << "\t\t******the above chain is not valid******" <<std::endl;
#endif

	return still_valid;
}

bool certificate::add_bg_path(const chain * c) throw() {
	list<node> l;
	chain_node_iterator it(c, decomposition);

	for(node n = it.next(); n!=NULL; n=it.next()) {
		l.append(n);
	}

	add_bg_path(l);
	return true;
}

bool certificate::verify() throw() {
	if (!still_valid)
		return false;
	if (!graphs_isomorphic(the_graph, new_graph, new_2_orig)) {
#ifndef NDEBUG
		std::cout << "resulting graphs not isomorphic";
#endif
		return false;
	}
	for(int i = chains.size()-1; i>=0; i--) {

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
#ifndef NDEBUG
				std::cout << "chain between parallel links"<< std::endl;
#endif
				return false; //chain between parallel links
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
#ifndef NDEBUG
				std::cout << "third node is one of the endpoints " << std::endl;
#endif
				return false; //third node is one of the endpoints
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

		{	edge to_be_deleted;
			forall(to_be_deleted, *chains[i]) {
				// delete the chain
				delete(le_edges[to_be_deleted]);
				le_edges[to_be_deleted] = NULL;

				new_graph.del_edge(to_be_deleted);
			}
			delete(chains[i]);
			chains[i] = NULL;
		}
	}

	return true;
}
