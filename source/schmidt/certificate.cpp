#include "certificate.hpp"
#include "chain_edge_iterator.hpp"
#include "utilities.hpp"

//#define BGCOUT
//#define VERIFYCOUT
//#define DO_NOTHING

certificate::certificate(ugraph  & graph,   schmidt_triconnectivity* d) :
	still_valid(true),
	the_graph(graph),
	decomposition(d),
	edge_items(new_graph,the_graph.number_of_edges()+the_graph.number_of_nodes(), NULL),
	orig_2_new(the_graph,NULL),
	edge_accounted_for(the_graph,false)
	{}

certificate::~certificate() {
	edge e;
	forall_edges(e,new_graph) {
		if (edge_items[e]!=NULL)
			delete(edge_items[e]);
	}
}

bool certificate::add_bg_path(list<edge> const & edges) throw() {
#ifdef DO_NOTHING
	return true;
#else
	if (!still_valid)
		return false;
#ifdef BGCOUT
	std::cout << "\t****** BG PATH: " ;
#endif
	node first_node = NULL;
	node last_node = NULL;
	node current_node = NULL;
	paths.push_back(new list<edge>());

	unsigned int prev_created_nodes = 0;

	node last=NULL,nodes[2]={NULL,NULL};
	unsigned int num=0;
	edge e;
	forall(e,edges) {
		edge_accounted_for[e]=true;
		if (last==NULL) {
			nodes[0] = source(e);
			nodes[1] = target(e);
			num = 2;

			last = target(e);
		} else {
			nodes[0] = opposite(last,e);
			nodes[1] = NULL;
			num=1;
			last = opposite(last,e);	//std::vector<std::pair<node,node>* > endvertices;
		}
		for(unsigned int i = 0; i<num; i++) {

			current_node = nodes[i];
			if (first_node == NULL)
				first_node = current_node;

			node new_node;
			if (orig_2_new[current_node] == NULL) {
				new_node = new_graph.new_node();
				orig_2_new[current_node] = new_node;
			} else {
				new_node = orig_2_new[current_node];
				prev_created_nodes++;
			}

			if (last_node !=NULL) {
				node prev_node = orig_2_new[last_node];
				assert(prev_node !=NULL);
				edge added_edge = new_graph.new_edge(prev_node,new_node);
				list<edge>::item it = paths.back()->append(added_edge);
				edge_items[added_edge] = new pair<list<edge>*,list<edge>::item>(paths.back(), it);
			}

			last_node = current_node;

		}
	}

	// paths have only two nodes in common with the rest of the graph
	still_valid &= paths.size() == 1 || prev_created_nodes == 2;

	// bg paths mustn't form cycles
	still_valid &= first_node != last_node;

#ifdef BGCOUT
	std::cout << std::endl;
#endif

#ifndef NDEBUG
	if (!still_valid)
		std::cout << "\t\t******the above chain is not valid******" <<std::endl;
#endif

	return still_valid;
#endif
}

bool certificate::add_bg_path(const chain * c) throw() {
#ifdef DO_NOTHING
	return true;
#else
	list<edge> l;
	chain_edge_iterator it(c, decomposition);

	for(edge n = it.next(); n!=NULL; n=it.next()) {
		l.append(n);
	}

	return add_bg_path(l);
#endif
}

bool certificate::verify() throw() {
#ifdef DO_NOTHING
	return true;
#else
	if (!still_valid)
		return false;
	bool isomorphic = true;
	{edge e;
	forall_edges(e,the_graph) {
		isomorphic &= edge_accounted_for[e];
	}}
	if (!isomorphic)
		return false;

	for(int the_path = paths.size()-1; the_path>=3; the_path--) {
#ifdef VERIFYCOUT
		std::cout << "Checking path " << the_path << std::endl;
		{	node n;
			forall_nodes(n,new_graph) {
				std::cout <<"(" << n->id() << "," << new_graph.degree(n) << ") ";
			}
			std::cout << std::endl;

		}
#endif
		list<edge>* current_path =paths[the_path];

		if (current_path->size()!=1)
			return false;

		const edge bg_edge = current_path->Pop();
		delete(current_path);
		paths[the_path]=NULL;

		const node a = source(bg_edge);
		const node b = target(bg_edge);

		assert(graph_of(a) == &new_graph);
		assert(graph_of(b) == &new_graph);
		const bool a_smooth = new_graph.degree(a) == 3;
		const bool b_smooth = new_graph.degree(b) == 3;
		unsigned int type = 0;
		if (a_smooth) type++;
		if (b_smooth) type++;

		delete(edge_items[bg_edge]);
		edge_items[bg_edge] = NULL;
		assert(bg_edge!=NULL);
		assert(graph_of(bg_edge) == &new_graph);
		new_graph.hide_edge(bg_edge);

		switch(type) {
		case 2: {		//subdivide 2, connect
#ifdef VERIFYCOUT
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
					b_neighbours[i++]  = opposite(e,b);
				}
			}


			if (a_neighbours[0] < a_neighbours[1]) {
				std::swap(a_neighbours[0], a_neighbours[1]);
			}
			if (b_neighbours[0] < b_neighbours[1]) {
				std::swap(b_neighbours[0], b_neighbours[1]);
			}
			if (the_path>3 && a_neighbours[0] == b_neighbours[0] && a_neighbours[1] == b_neighbours[1]) {
#ifndef NDEBUG
				std::cout << "chain between parallel links"<< std::endl;
#endif
				return false;
			} else if(the_path==3 && !(a_neighbours[0] == b_neighbours[0] && a_neighbours[1] == b_neighbours[1])) {
#ifndef NDEBUG
				std::cout << "chain 3 does not run between parallel links" << std::endl;
#endif
				return false;
			}
			// update edge_items and the paths
#ifdef VERIFYCOUT
			std::cout << "updating" << std::endl;
#endif
			{ 	list<edge>* lists[2]={NULL,NULL};
				edge es[2] = {NULL,NULL};

				//first remove the smoothened edges from the paths
				{	unsigned int j=0;
					edge e;
					forall_adj_edges(e,a) {
						pair<list<edge>*, list<edge>::item>* p = edge_items[e];
						p->first->del_item(p->second);
						lists[j] = p->first;
						delete(edge_items[e]);
						edge_items[e]=NULL;
					}
					es[0] = smoothen(new_graph,a);
					j++;
					forall_adj_edges(e,b) {
						pair<list<edge>*, list<edge>::item>* p = edge_items[e];
						p->first->del_item(p->second);
						lists[j] = p->first;
						delete(edge_items[e]);
						edge_items[e]=NULL;
					}
				}

				es[1] = smoothen(new_graph,b);

				//and insert the new edges into the paths
				for(unsigned int i = 0; i<2; i++) {
					list<edge>::item it = lists[i]->append(es[i]);
					assert(es[i] != NULL);
					assert(graph_of(es[i]) == &new_graph);
					edge_items[es[i]] = new pair<list<edge>*,list<edge>::item>(lists[i],it);
				}

			}

		} break;
		case 1: {		//subdivide 1, connect
#ifdef VERIFYCOUT
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

			//update edge_items and the paths
			//first delete smoothened edges from the chain
#ifdef VERIFYCOUT
			std::cout << "updating" << std::endl;
#endif

			{ 	list<edge>* l=NULL;

				{	edge e;
					forall_adj_edges(e,sub) {
						pair<list<edge>*, list<edge>::item>* p = edge_items[e];
						p->first->del_item(p->second);
						l = p->first;
						delete(edge_items[e]);
						edge_items[e]=NULL;
					}
				}

				//then from the graph
				edge e1 = smoothen(new_graph,sub);

				//and insert the new edge into the chain
				list<edge>::item it = l->append(e1);

				assert(e1!=NULL);
				assert(graph_of(e1) == &new_graph);

				edge_items[e1] = new pair<list<edge>*, list<edge>::item>(l,it);
			}
		} break;
		case 0: {		//connect two is always okay, because of the checks we performed before adding this chain.
#ifdef VERIFYCOUT
			std::cout << "connect two" << std::endl;
#endif
		} break;
		default:
			assert(false);
		}
	}


	node remaining_nodes[2]={NULL,NULL};
	for(int the_path=2; the_path>=0; the_path-- ) {
		list<edge>* current_path = paths[the_path];
		edge bg_edge = current_path->Pop();
		delete(paths[the_path]);
		paths[the_path]=NULL;
		delete(edge_items[bg_edge]);
		edge_items[bg_edge] = NULL;
		node endvertices[2];
		endvertices[0] = source(bg_edge);
		endvertices[1] = target(bg_edge);
		if (endvertices[0] < endvertices[1])
			std::swap(endvertices[0], endvertices[1]);
		for(unsigned int i = 0; i<2; i++)
			if (remaining_nodes[i] == NULL)
				remaining_nodes[i] = endvertices[i];

		for(unsigned int i=0; i<2; i++)
			if (endvertices[i] != remaining_nodes[i])
				return false;

		new_graph.hide_edge(bg_edge);
	}
	return true;
#endif
}
