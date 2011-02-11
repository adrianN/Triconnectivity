#include "certificate.hpp"
#include "chain_edge_iterator.hpp"
#include "utilities.hpp"

#define BGCOUT

certificate::certificate(ugraph  & graph,   schmidt_triconnectivity* d) :
	still_valid(true),
	the_graph(graph),
	decomposition(d),
	created_by_chain(new_graph,the_graph.number_of_nodes(), -1),
	le_edges(new_graph,the_graph.number_of_edges(), NULL),
	orig_2_new(the_graph,NULL),
	new_2_orig(new_graph,the_graph.number_of_nodes(),NULL),
	edge_accounted_for(the_graph,false)
	{}

certificate::~certificate() {
}

bool certificate::add_bg_path(list<edge> const & edges) throw() {
	if (!still_valid)
		return false;
#ifdef BGCOUT
	std::cout << "\t****** BG PATH: " ;
#endif
	node first_node = NULL;
	node last_node = NULL;
	node current_node = NULL;
	chains.push_back(new list<edge>());

	unsigned int prev_created_nodes = 0;

	node last=NULL,nodes[2]={NULL,NULL};
	unsigned int num=0;
	list<node> l;
	edge e;
	forall(e,edges) {
		edge_accounted_for[e]=true;
		if (last==NULL) {
			nodes[0] = source(e);
			nodes[1] = target(e);
			num = 2;
			l.append(source(e));
			l.append(target(e));
			last = target(e);
		} else {
			l.append(opposite(last,e));
			nodes[0] = opposite(last,e);
			nodes[1] = NULL;
			num=1;
			last = opposite(last,e);
		}
		for(unsigned int i = 0; i<num; i++) {

			current_node = nodes[i];
	#ifdef BGCOUT
			std::cout << decomposition->dfi(current_node) << " ";
	#endif
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


	}

	// paths have only two nodes in common with the rest of the graph
	still_valid &= chains.size() == 1 || prev_created_nodes == 2;

	// bg paths mustn't form cycles
	still_valid &= first_node != last_node;

	endvertices.push_back(new pair<node,node>(orig_2_new[first_node],orig_2_new[last_node]));
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
	list<edge> l;
	chain_edge_iterator it(c, decomposition);

	for(edge n = it.next(); n!=NULL; n=it.next()) {
		l.append(n);
	}

	add_bg_path(l);
	return true;
}

bool certificate::verify() throw() {
	if (!still_valid)
		return false;
	bool isomorphic = true;
	{edge e;
	forall_edges(e,the_graph) {
		isomorphic &= edge_accounted_for[e];
	}}
	if (!isomorphic)
		return false;
	if (!graphs_isomorphic(the_graph, new_graph, new_2_orig)) {
#ifndef NDEBUG
		std::cout << "resulting graphs not isomorphic";
#endif
		return false;
	}
	for(int i = chains.size()-1; i>=3; i--) {
		{	node n;
		forall_nodes(n,new_graph) {
			std::cout <<"(" << n->id() << "," << new_graph.degree(n) << ") ";
		}
		std::cout << std::endl;

	}
		const node a = endvertices[i]->first;
		const node b = endvertices[i]->second;
		delete endvertices[i];
		endvertices[i] = NULL;

		assert(graph_of(a) == &new_graph);
		assert(graph_of(b) == &new_graph);
		const bool a_smooth = new_graph.degree(a) == 3;
		const bool b_smooth = new_graph.degree(b) == 3;
		unsigned int type = 0;
		if (a_smooth) type++;
		if (b_smooth) type++;
#ifdef BGCOUT
		std::cout << "Checking chain " << i << " endvertices " << decomposition->dfi(new_2_orig[a]) << " (" << new_graph.degree(a) << "," << a->id() << ") " << decomposition->dfi(new_2_orig[b]) << " (" << new_graph.degree(b) << "," << b->id() << ") " << type<< " ";
		std::cout << "edges ";
		{ edge e;
			forall(e, *chains[i]) {
				std::cout << decomposition->dfi(new_2_orig[source(e)]) << " " << decomposition->dfi(new_2_orig[target(e)]) << " " ;
			}
			std::cout << std::endl;
		}
#endif

		{	edge to_be_deleted;
			forall(to_be_deleted, *chains[i]) {
				// delete the chain
				delete(le_edges[to_be_deleted]);
				le_edges[to_be_deleted] = NULL;

				new_graph.del_edge(to_be_deleted);
				if (new_graph.degree(source(to_be_deleted)) == 0) {
					assert(created_by_chain[source(to_be_deleted)] == i);
				}
				if (new_graph.degree(target(to_be_deleted)) == 0) {
					assert(created_by_chain[target(to_be_deleted)] == i);
				}

			}
			delete(chains[i]);
			chains[i] = NULL;
		}

		switch(type) {
		case 2: {		//subdivide 2, connect
#ifdef BGCOUT
			std::cout << "subdivide 2, connect" << std::endl;
#endif
			node a_neighbours[2];
			node b_neighbours[2];
			{	edge e;
				unsigned int i = 0;
				forall_adj_edges(e,a) {
					assert(graph_of(e) == &new_graph);
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
			if (i!=3 && a_neighbours[0] == b_neighbours[0] && a_neighbours[1] == b_neighbours[1]) {
#ifndef NDEBUG
				std::cout << "chain between parallel links"<< std::endl;
#endif
				return false;
			} else if(i==3 && !(a_neighbours[0] == b_neighbours[0] && a_neighbours[1] == b_neighbours[1])) {
				std::cout << "chain 3 does not run between parallel links" << std::endl;
				return false;
			}
			// update le_edges and the chains
			std::cout << "updating" << std::endl;
			{ 	list<edge>* lists[2][2];

				//first remove the smoothened edges from the chains
				{	unsigned int i=0, j=0;
					edge e;
					forall_adj_edges(e,a) {
						pair<list<edge>*, list<edge>::item>* p = le_edges[e];
						p->first->del_item(p->second);
						lists[j][i++] = p->first;
					}
					assert(i==2);
					j++;
					i=0;
					forall_adj_edges(e,b) {
						pair<list<edge>*, list<edge>::item>* p = le_edges[e];
						p->first->del_item(p->second);
						lists[j][i++] = p->first;
					}
					assert(j==1);
				}

				// then remove them from the graph
				edge es[2];
				es[0] = smoothe(new_graph,a);
				es[1] = smoothe(new_graph,b);

				//and insert the new edges into the chains
				for(unsigned int i = 0; i<2; i++) {
					std::cout << i << " belonged to " << lists[i][0] << " " << lists[i][1] << std::endl;
					assert(lists[i][1] == lists[i][0]);
					list<edge>::item it = lists[i][0]->append(es[i]);
					le_edges[es[i]] = new pair<list<edge>*,list<edge>::item>(lists[i][0],it);
				}
			}

		} break;
		case 1: {		//subdivide 1, connect
#ifdef BGCOUT
			std::cout << "subdivide one, connect" << std::endl;
#endif
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
			std::cout << "updating" << std::endl;

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
#ifdef BGCOUT
			std::cout << "connect two" << std::endl;
#endif
		} break;
		default:
			assert(false);
		}
	}

	//TODO check 0,1,2

	return true;
}
