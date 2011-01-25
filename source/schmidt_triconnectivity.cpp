#include "schmidt.hpp"
#include "dfs.hpp"
#include "InterlockingIntervals.hpp"
#include "chain_node_iterator.hpp"
#include "LEDA/core/stack.h"
#include "LEDA/graph/ugraph.h"

using namespace leda;

#define RECORD_NODES

schmidt_triconnectivity::schmidt_triconnectivity(ugraph& graph) :
	belongs_to_chain(graph,-1),
	is_frond(graph, false),
	dfis(graph,0),
	node_at(new node[graph.number_of_nodes()+1]),
	parent(graph, NULL),
	inner_of_chain(graph,-1),
	type_3(graph),
	is_real(graph,false),
	chains(new chain[graph.number_of_edges() - graph.number_of_nodes() + 2]),
	caterpillars(new caterpillar[graph.number_of_edges() - graph.number_of_nodes() +2]), //TODO not space efficient, at most #fronds/2 caterpillars
	the_graph(graph)
{
	{	edge e;
		forall_edges(e,the_graph) {
			belongs_to_chain[e] = -1; //leda bugfix
		}
	}
	{ 	node n; //leda bugfix
		forall_nodes(n,the_graph) {
			inner_of_chain[n] = -1;
			is_real[n] = false;
		}
	}

#ifndef NDEBUG
	for(unsigned int i =0; i<(unsigned int) the_graph.number_of_nodes(); i++) node_at[i]=NULL;
#endif
	assert(is_connected(the_graph));
	assert(the_graph.number_of_nodes()>=4);
	initial_dfs();
	chain_decomposition();

	//Assert property A TODO handle without asserts
	const unsigned int number_chains = graph.number_of_edges() - graph.number_of_nodes() + 1;
	node min_degree_vertex = the_graph.first_node();
	{	node n;
		forall_nodes(n,the_graph) {
			if (the_graph.degree(n) < the_graph.degree(min_degree_vertex))
				min_degree_vertex = n;
		}
	}
	switch (the_graph.degree(min_degree_vertex)) {
	case 1: assert(false && "not biconnected"); break;
	case 2: assert(false && "not triconnected"); break;
	case 3: {
		//count cycles in the chain decomposition
		unsigned int num_cycles = 0;
		for(unsigned int i = 0; i<=number_chains; i++) {
			if (chains[i].s == chains[i].t) num_cycles++;
		}
		assert(num_cycles==0); //should already happen during the decomposition. Not 1 as in thesis, cause C0
		break;
	}
	default: break;
	}

}

schmidt_triconnectivity::~schmidt_triconnectivity(void) {
	delete[] chains;
	delete[] node_at;
	delete[] caterpillars;
}

certificate* schmidt_triconnectivity::certify(void) {
	for(unsigned int i=1; i<(unsigned int)(the_graph.number_of_edges() - the_graph.number_of_nodes() + 2); i++) {
		const chain& current_chain = chains[i];
		assert(in_subdivision(current_chain));

		//Compute children 12, the set of all children of type 1 or 1 of current_chain that are not yet added to the subdivision
		slist<chain*> children12;
		int_set set_children12((the_graph.number_of_edges() - the_graph.number_of_nodes() + 2));
		{	chain* child;
			forall(child,children12) {
				if (!in_subdivision(child) && (child->type == one || child->type == two_a || child->type == two_b)) {
					children12.append(child);
					set_children12.insert(child->number);
				}
			}
		}

		//compute type3, the set of chains of type 3 that start at an inner vertex of current_chain
		slist<chain*> type3;
		for(node v = current_chain.t; dfi(v) < dfi(current_chain.s); v = parent_node(v)) { //<=?
			type3.conc(*type_3[v]); //TODO order?
		}

		/* preemptively add chains of type 2a with real t(C). As they are backedges with parent in S they have prop. 30a. As t(C)
		 * is real and a proper descendant of t(current_chain) (which is real) by definition of type 2 chains, they have property 30b.
		 */
		preemptively_add_easy_clusters(children12);


		//Partition type3 into segments. The key is the number of the minimal chain in the segment
		h_array<unsigned int, slist<chain*> > segments;
		// Also compute the attachment vertices of each segment, although they aren't needed if it is an easy cluster
		// The attachment vertices are the end vertices of the minimal chain and of the chains in type3 that are in the segment
		h_array<unsigned int, slist<node> > attachment_vertices;
		//true if a segment has no chain in children12
		h_array<unsigned int, bool> intersection_12_free(true);

		partition_into_segments(current_chain, type3, segments, attachment_vertices, set_children12, intersection_12_free);

		//Add the easy segments, i.e. those without chains from children12
		add_easy_segments(current_chain, type3,intersection_12_free);

		// Add the hard clusters
		// we take the attachment vertices and generate a set of intervals for each segment. Then we compute an order
		// for the intervals such that we can add the minimal chain (and thus all others in the segment) such that there is at least one real vertex between
		// the end nodes of the chain.

		//add_hard_segments();


	}
	return NULL;
}

void schmidt_triconnectivity::preemptively_add_easy_clusters(slist<chain*> children12) {
	{
		slist<chain*>::item prev_it=NULL;
		for(slist<chain*>::item cur_it = children12.first(); cur_it!=NULL; (prev_it=cur_it, cur_it=children12.next_item(cur_it))){
			chain* cur_child = children12.contents(cur_it);
			if (cur_child->type == two_a && is_real[cur_child->t]) {
				add_to_subdivision(cur_child);
				//delete from list
				if (prev_it!=NULL) {
					children12.del_succ_item(prev_it);
					cur_it=prev_it;
				} else {
					children12.pop();
					//TODO does next_item work as expected? Otherwise switch to doubly linked lists. Otherwise another pass over the list...
				}
			}
		}
	}
}

void schmidt_triconnectivity::partition_into_segments(const chain& current_chain, const slist<chain*>& type3, h_array<unsigned int, slist<chain*> >& segments, h_array<unsigned int, slist<node> >& attachment_vertices, const int_set& set_children12, h_array<unsigned int,bool>& intersection_12_free) {
	{	node_array<int> minimal_chain(the_graph,-1);
		{	node n; //leda bugfix
			forall_nodes(n,the_graph) {
				minimal_chain[n] = -1;
			}
		}
		chain* t3_chain;
		forall(t3_chain, type3) {
			//Find the minimal chain
			node v = t3_chain->t;
			do {
				v = parent_node(v);
			} while (!in_subdivision(parent_node(v)) && minimal_chain[v] < 0);

			if (minimal_chain[v] < 0) {
				minimal_chain[v] = inner_of_chain[parent_node(v)];
			}

			assert(minimal_chain[v] >= 0);
			unsigned int m_chain = (unsigned int)minimal_chain[v];

			assert(contained_in_chain(chains[m_chain].s, &current_chain));
			assert(contained_in_chain(chains[m_chain].t, &current_chain));


			//mark all nodes on the path with the minimal chain to get linear running time. Unsure whether we have
			//to keep the markings between different C_i. Probably not.
			v = t3_chain->t;
			do {
				minimal_chain[v] = m_chain;
				v = parent_node(v);
			} while (!in_subdivision(parent_node(v)) && minimal_chain[v] < 0);

			t3_chain->segment = m_chain;
			assert(segments[m_chain].contents(segments[m_chain].last())->number < t3_chain->number);
			if (!segments.defined(m_chain)) {
				//we add a new segment, put the end vertices of the minimal chain in the set of attachment vertices
				chain& min_chain = chains[m_chain];
				attachment_vertices[m_chain].append(min_chain.s);
				attachment_vertices[m_chain].append(min_chain.t);
			}
			segments[m_chain].append(t3_chain);
			attachment_vertices[m_chain].append(t3_chain->s);
			attachment_vertices[m_chain].append(t3_chain->t);

			intersection_12_free[m_chain] &= !set_children12.member(m_chain);
		}
	}
}

void schmidt_triconnectivity::add_easy_segments(const chain& current_chain, slist<chain*>& type3, const h_array<unsigned int, bool>& intersection_12_free) {

	{ 	slist<chain*>::item prev_it=NULL;
		for(slist<chain*>::item cur_it = type3.first(); cur_it!=NULL; (prev_it=cur_it, cur_it=type3.next_item(cur_it))){
			chain* t3_chain = type3.contents(cur_it);
			assert(t3_chain->segment>=0);

			if (!intersection_12_free[t3_chain->segment]) continue;

			//add all ancestors of t3_chain that are still in H, in order of <
			stack<chain*> ancestors_in_segment;
			for(chain* cur_chain=t3_chain; cur_chain->parent != NULL && cur_chain->number > (unsigned int)t3_chain->segment; cur_chain = cur_chain->parent) {
				//make sure cur_chain is not a type 3 chain that won't be removed from type3(C_i). I hope this never happens
				assert(cur_chain==t3_chain || cur_chain->type != three_a || cur_chain->type != three_b || !contained_in_chain(cur_chain->s, &current_chain));

				ancestors_in_segment.push(cur_chain);
			}

			while(ancestors_in_segment.size()>0) {
				add_to_subdivision(ancestors_in_segment.pop());
			}

			if (prev_it!=NULL) {
				type3.del_succ_item(prev_it);
				cur_it=prev_it;
			} else {
				type3.pop();
				//TODO does next_item work as expected?
			}
		}
	}
}

void schmidt_triconnectivity::add_hard_segments(const chain& current_chain, h_array<unsigned int, slist<node> > attachment_vertices) {
	//we start by computing an order of the segments
	// we map the nodes of the current_chain to ints from 1 to |current_chain|
	// we add an artificial interval for each real node on the chain

	slist<interval<chain* > * > intervals;
	node_array<unsigned int> mapping(the_graph);
	std::vector<node> reverse_mapping;
	reverse_mapping.push_back(NULL); //mapping starts at 1
	unsigned int vertices_in_chain = 0;

	chain_node_iterator it(current_chain, this);
	for(node n=it.next(); it.has_next(); n=it.next()) {
		vertices_in_chain++;
		mapping[n] = vertices_in_chain;
		reverse_mapping.push_back(n);

		//create artificial intervals for vertices that are already real
		if (is_real[n]) {
			interval<chain*>* val  = new interval<chain*>(0,vertices_in_chain, NULL); // artificial intervals for real vertices start at 0
			intervals.append(val);
		}
	}

	//iterate over the attachment vertices of all chains and generate the intervals
	//we generate one interval from the lowest attachment vertex to all others and one from the highest to all others
	slist<slist<interval<chain*>* >* > equivalent_intervals;
	unsigned int m_chain;
	forall_defined(m_chain, attachment_vertices) {
		slist<node>& vertices = attachment_vertices[m_chain];
		slist<unsigned int> mapped_vertices;
		node n;
		//TODO this can be simplified by walking in the right direction and avoiding one bucket sort
		forall(n,vertices) {
			mapped_vertices.append(mapping[n]);
		}

		{
			unsigned int (*id)(const unsigned int& a) = (identity<unsigned int>);
			mapped_vertices = bucket_sort<unsigned int>(mapped_vertices, id, 0, vertices_in_chain);
		}

		slist<interval<chain*>*>* equivalence_class = new slist<interval<chain*>*>();
		equivalent_intervals.append(equivalence_class);

		//generate the intervals
		const unsigned int smallest = mapped_vertices.contents(mapped_vertices.first());
		const unsigned int largest = mapped_vertices.contents(mapped_vertices.last());
		interval<chain*>* val = new interval<chain*>(smallest,largest,&chains[m_chain]);
		intervals.append(val);
		equivalence_class->append(val);

		slist<unsigned int>::item it = mapped_vertices.first();
		it = mapped_vertices.succ(it);

		for(unsigned int i=1; i <= (unsigned int)mapped_vertices.size()-2; i++) {
			val = new interval<chain*>(smallest, mapped_vertices.contents(it), &chains[m_chain]);
			intervals.append(val);

			val = new interval<chain*>(mapped_vertices.contents(it), largest, &chains[m_chain]);
			intervals.append(val);
			equivalence_class->append(val);
		}
	}

	// sort the intervals

	unsigned int (*fc)(interval<chain*>* const &) = &interval<chain*>::first_component;
	slist<interval<chain*>*> intervals_asc = bucket_sort<interval<chain*>*>(intervals,fc,0,vertices_in_chain);
	unsigned int (*sc)(interval<chain*>* const &) = &interval<chain*>::second_component;
	slist<interval<chain*>*> intervals_dsc = bucket_sort<interval<chain*>*>(intervals,sc,0,vertices_in_chain);


	// order the intervals in fancy order



	// free resources
	{
		interval<chain*>* val;
		forall(val, intervals_asc) {
			delete val;
		}
		slist<interval<chain*>*>* equivalence_class;
		forall(equivalence_class, equivalent_intervals) {
			delete equivalence_class;
		}
	}

	// add the clusters in the order of the intervals

}

inline unsigned int schmidt_triconnectivity::dfi(const node v) const {return dfis[v];}
inline void schmidt_triconnectivity::set_dfi(const node v, const unsigned int d) { assert(dfis[v] == 0); dfis[v] = d; node_at[d] = v;}

node schmidt_triconnectivity::parent_node(const node v) const {
	edge e = parent[v];
	assert(!is_frond[e]);
	assert(source(e)==v);
	return target(e);
}

inline bool schmidt_triconnectivity::in_subdivision(const chain& c) const {
	return c.in_subdivision;
}

inline bool schmidt_triconnectivity::in_subdivision(const chain* c) const {
	return c->in_subdivision;
}

inline bool schmidt_triconnectivity::in_subdivision(const node n) const {
	return n == node_at[1] || inner_of_chain[n] >= 0; //every node is the inner node of some chain, except maybe the root. I'm not sure about that TODO
}

void schmidt_triconnectivity::chain_decomposition(void) {
	//Calculate the list of fronds in DF order of their starting vertex and find the chain for each frond. This can probably be incorporated into the DFS.
	unsigned int chain_number = 3; //the first three chains are found during dfs

	for(unsigned int i=1; i<(unsigned int)the_graph.number_of_nodes()+1; i++) {
		node current_node = node_at[i];
		assert(current_node!=NULL);
		edge e;
		forall_adj_edges(e,current_node) {
			if (!is_frond[e] || belongs_to_chain[e] >= 0) // the first three chains already exist
				continue;

			chain& current_chain = chains[chain_number];

			//fronds go from up to down, hence chain.s = source, chain.t = target
			belongs_to_chain[e] = chain_number;
			assert(current_chain.s == NULL);
			current_chain.s = source(e);

			mark_path(target(e),chain_number); //sets .t
			chains[chain_number].first_edge = e;

			assert(current_chain.s!=current_chain.t && "graph is not biconnected");
			assert(dfi(current_chain.s) < dfi(current_chain.t));

			switch (classify_chain(chain_number)) {
			case three_a: case three_b:
				type_3[current_chain.s]->append(&current_chain);
			default: break;
			}

			chain_number++;
		}
	}



	std::fstream f("./blubb.dot", std::ios::out);
	dfs_tree_to_dot(f);
	f.close();

#ifndef NDEBUG
	//Chain decompositition partitions edges
	{	edge e;
		forall_edges(e,the_graph) {
			assert(belongs_to_chain[e] >= 0);
		}
	}
#endif

}

/* In this function we do a DFS on the graph to
 * * calculate DFIs for each vertex
 * * finding an initial K32 subdivision
 */
void schmidt_triconnectivity::initial_dfs(void) {


		node current_node = source(the_graph.first_edge());
		const node root = current_node;
		inner_of_chain[root] = 0;
		assert(the_graph.degree(root)>=3);

		unsigned int number_seen=1;
		set_dfi(current_node, number_seen++);

		parent[current_node] = NULL;
		edge_array<bool> seen_edge(the_graph, false);

		chains[0].parent = NULL;
		chains[0].number = 0;
		chains[1].parent = &chains[0];
		chains[2].parent = &chains[0];
		chains[1].s = root;
		chains[2].s = root;

		/* first find the first fundamental cycle back to the root
		 * do a normal dfs
		 */
		std::cout << "first cycle " << std::endl;
		if (!find_a_cycle_to_root(current_node, number_seen, seen_edge, false,  1)) {
			assert(false && "graph is not triconnected");
		}

		const node node_a = current_node;

		std::cout << "node a: "  << current_node->id() << std::endl;
		/* When the first cycle is found we start a new DFS for every node on the path from a upwards to the root. When one
		 * of those searches encounters another backedge to the root we know that the starting node is the LCA.
		 *
		 * TODO we can also start labeling edges as belonging to a path
		 */

		node inner_search_cur_node;
		node lca = NULL;
		node node_b = NULL;
		while(current_node != NULL) { //todo switch from node to edge
			std::cout << "walk up " << current_node->id() <<std::endl;

			inner_search_cur_node = current_node;

			if (find_a_cycle_to_root(inner_search_cur_node,number_seen, seen_edge, true, 2) && lca == NULL && node_b == NULL) {
				std::cout << "found a cycle" << std::endl;
				lca = current_node;
				node_b = inner_search_cur_node;

			}
			if (parent[current_node]!=NULL) {
				std::cout << "the parent of " << dfi(current_node) << " is " << parent[current_node] << " " << dfi(source(parent[current_node])) << " " << dfi(target(parent[current_node])) << std::endl;
				assert(source(parent[current_node]) == current_node);
				current_node = target(parent[current_node]);
			} else {
				current_node = NULL;
			}

		}
		assert(lca!=NULL);
		assert(node_b!=NULL);
		std::cout << "Root " << root->id() << std::endl;
		std::cout << "LCA " << lca->id() << std::endl;

#ifndef NDEBUG
		//DFS visits all edges
		{	edge e;
			forall_edges(e,the_graph) {
				assert(seen_edge[e] && "unvisited edge");
				if (!((is_frond[e] && dfi(source(e))<dfi(target(e))) || (!is_frond[e] && dfi(source(e)) > dfi(target(e))))) {
					std::cout << "edge " << e << " " << dfi(source(e)) << " " << dfi(target(e)) << " is frond? " << is_frond[e] << std::endl;
					assert(false);
				}
			}
		}
#endif

		//TODO, do we need these at all? C0 stays unclassified
		chains[0].s = lca;
		mark_path(lca, 0);

		std::cout << "A " << node_a->id() << std::endl;
		add_to_subdivision(&chains[0]);


		mark_path(node_a, 1);
		std::cout << "B " << node_b->id() << std::endl;
		add_to_subdivision(&chains[1]);

		mark_path(node_b, 2);
		add_to_subdivision(&chains[2]);
		assert(chains[0].s == chains[2].t);


		for(unsigned int i= 0 ; i<3; i++)
			classify_chain(i);



}

/* takes the first inner vertex of a chain (except C0?) and walks upwards in the tree until the edge to the parent is contained in another chain
 * Marks edges on the path as belonging to the current chain. Marks inner nodes as inner nodes of the chain. Sets chain.t */
void schmidt_triconnectivity::mark_path(const node start, const unsigned int the_chain) {
	chains[the_chain].number = the_chain;

	node current_node = start;
	for(edge current_edge = parent[start]; current_edge!=NULL; (current_node = parent_node(current_node), current_edge = parent[current_node])) {
		chains[the_chain].t = current_node;

		if (belongs_to_chain[current_edge]>=0) {
			//found an edge that belongs to a different chain, hence this chain ends here.
			assert((unsigned int)belongs_to_chain[current_edge] < the_chain && "the parent of a chain must be a chain of smaller chain number");
			chains[the_chain].parent = &chains[belongs_to_chain[parent[current_node]]];
			break;
		}

		chains[the_chain].is_backedge = false; // more than one edge contained

		assert(inner_of_chain[current_node] < 0 && "No node is inner node of >1 chains");
		inner_of_chain[current_node] = the_chain;
		//must be <0 due to if
		belongs_to_chain[current_edge] = the_chain;
	}
}

inline bool schmidt_triconnectivity::contained_in_chain(const node v, const chain* chain) {
	return (unsigned int)inner_of_chain[v] == chain->number || v == chain->s || v == chain->t;
}

chain_type schmidt_triconnectivity::classify_chain(const unsigned int chain_number) {
	if (chain_number == 0)
		return unmarked; // first chain stays unmarked

	chain* const the_chain = &chains[chain_number];
	chain* const the_parent = the_chain->parent;

	assert(contained_in_chain(the_chain->t,the_parent));

	// t(the_chain) ->T s(the_chain) is contained in parent
	if (the_parent->number == 0 || (the_chain->s != the_parent->s && contained_in_chain(the_chain->s, the_parent))) {
#ifndef NDEBUG
		//walk the path the_chain.t ->T the_chain.s and check if each node is inner_node in parent
		{
			node current_node = the_chain->t;
			assert((unsigned int)inner_of_chain[current_node] == the_parent->number);
			do {
				current_node = parent_node(current_node);
				assert(contained_in_chain(current_node,the_parent));
			} while (current_node != the_chain->s);
		}
#endif
		the_chain->type = one;
		return one;
	}

	//the chain starts at the same node as its parent, it is of type 2
	if (the_chain->s == the_parent->s) {
		assert(the_parent->number != 0);

		// type 2: the_chain = C_0, t(the_chain) is inner vertex of parent
		if (the_chain->is_backedge) {
			the_chain->type = two_a;
			return two_a;
		} else {
			the_chain->type = two_b;
			the_chain->is_marked = true;
			return two_b;
		}

	}


	if (!the_parent->is_marked) {
		the_chain->type = three_a;
		return three_a;
	} else {
		the_chain->type = three_b;
		caterpillars[chain_number].append(the_chain);
		chain* c_jay = the_parent;
		while(c_jay->is_marked) {
			assert(c_jay->type = two_b);

			c_jay->is_marked = false;
			caterpillars[chain_number].append(c_jay);
			c_jay = c_jay->parent;
		}
		return three_b;
	}
}

bool schmidt_triconnectivity::find_a_cycle_to_root(node& start_node, unsigned int& number_seen, edge_array<bool>& seen_edge, bool continue_after_found, const unsigned int current_chain) {
		stack<edge> first_cycle;

		node current_node = start_node;
		{	edge e;
			forall_adj_edges(e,current_node) {
				if (seen_edge[e])
					continue;

				first_cycle.push(e);
			}
		}

		bool found_cycle = false;
		while(!first_cycle.empty()) {
			const edge current_edge = first_cycle.pop();

			if (seen_edge[current_edge]) {
				continue;
			}

			seen_edge[current_edge] = true;

			//make sure fronds go from up to down and tree edges from down to up
			current_node = target(current_edge);
			node next_node = source(current_edge);

			assert(dfi(next_node) != 0 || dfi(current_node)!=0); //one of them must have been visited already

			//we don't know which way round the edges go :/. The next node must always have a lower dfi:
			//either it is unvisited, then it has dfi 0, or we have a frond, then dfi[next_node] < dfi[current_node] and both !=0

			if (dfi(next_node) > dfi(current_node)) {
				std::swap(next_node, current_node);
				the_graph.rev_edge(current_edge);
			}
			assert(dfi(current_node)!=0);
			assert(dfi(source(current_edge)) < dfi(target(current_edge)));
			assert(source(current_edge) == next_node);

			// mark fronds
			if (dfi(next_node)>0) {

				is_frond[current_edge] = true;

				if (dfi(next_node) == 1 && !found_cycle) { //backedge to the root
					//first cycle complete
					belongs_to_chain[current_edge] = current_chain;
					start_node = current_node; //pass current node outside
					found_cycle=true;
					if (!continue_after_found) {
						return found_cycle;
					}
				}
				continue;
			}
			//continue exploration
#ifndef NDEBUG
			if (parent[next_node]!=NULL) {
				std::cout << "Node " << next_node->id() << " already has parent " << source(parent[next_node])->id() << " " << target(parent[next_node])->id() << std::endl;
				assert(false);
			}
#endif

			parent[next_node] = current_edge;
			current_node = next_node;

			set_dfi(current_node,number_seen++);
			assert(dfi(source(current_edge)) > dfi(target(current_edge)));

			{	edge next_edge;
				//push unvisited neighbours
#ifndef NDEBUG
				int deg = 0;
#endif
				forall_adj_edges(next_edge,current_node) {
					assert(++deg <= the_graph.degree(current_node));

					//skip edge to parent
					if (seen_edge[next_edge]) {
						continue;
					}

					first_cycle.push(next_edge);
				}
				assert(deg == the_graph.degree(current_node));
			}
		}

		return found_cycle;
}

void schmidt_triconnectivity::add_to_subdivision(chain* c) {
	assert(!c->in_subdivision);
	assert((c->parent == NULL && c->number==0) || c->parent->in_subdivision); //modularity
	c->in_subdivision = true;
	is_real[c->s] = true;
	is_real[c->t] = true;

#ifndef NDEBUG
	if (c->number > 1) {
	edge e;
	unsigned int count=0;
		forall_adj_edges(e,c->s) {
			if (chains[belongs_to_chain[e]].in_subdivision) {
				count++;
			}
		}
		assert(count>=3);

		count=0;
		forall_adj_edges(e,c->t) {
			if (chains[belongs_to_chain[e]].in_subdivision) {
				count++;
			}
		}
		assert(count>=3);
	}
#endif
}

void schmidt_triconnectivity::dfs_tree_to_dot(std::ostream& out) {
	out << "digraph G {" << std:: endl;
	{
		node n;
		forall_nodes(n,the_graph) {
			unsigned int num = dfi(n);
			int id = n->id();
			/*
			 * +--------------------------------+
			 * | Parent             			|
			 * +--------------------------------+
			 * | LP1 | LP2 | highpoint | #desc. |
			 * +--------------------------------+
			 * | Node               			|
			 * +--------------------------------+
			 */
			out << "subgraph cluster" << inner_of_chain[n] << " { node" << id << "\n label = \"chain " << inner_of_chain[n] << "\" }" << std::endl;
			out << "\tnode";
			out << id;

#ifdef RECORD_NODES
			out << " [shape = record, label = \"{";
			if (parent[n] == NULL) {
				out << -1 << "|";
			} else {
				edge e = parent[n];

				out  << source(e)->id() << "|";
			}

			out << "Node " << num << ", ID " << id << " | ";
			out << inner_of_chain[n] << "}\"]" << std::endl;
#else
			out << " [label = \"";
			if (parent[n] == NULL) {
				out << -1 << "\\n";
			} else {
				edge e = parent[n];

				out  << source(e)->id() << "\\n";
			}
			out << num << " " << id << "\\n" << inner_of_chain[n] << \"]" << std::endl;
#endif
		}
	}
	{
		/* Edges are labeled with 1 if they start a path (0 oth.) and the value they by which they get ordered (phi) */
		edge e;
		forall_edges(e,the_graph) {
			node n = target(e);
			node u = source(e);

			out << "\tnode" << n->id() << " -> " << "node" << u->id();
			out << " [ " << (is_frond[e] ? "constraint = false, ":"") << "label=\"" << belongs_to_chain[e] <<  "\"]";
			out << std::endl;
		}

	}


	out << "}" << std::endl;
}

void schmidt_triconnectivity::chain_tree_to_dot(std::ostream& out) {
   out << "graph G {" << std::endl;
   const unsigned int num_chains = (unsigned int)(the_graph.number_of_edges() - the_graph.number_of_nodes()) +2;
   for(unsigned int i = 0; i<num_chains; i++) {
#ifdef RECORD_NODES
	   out << "node" << i << " [shape = record, label =\"{";
	   out << i << "| {";
	   out << "s " << (chains[i].s ? dfi(chains[i].s) : -1) << "| ";
	   out << "t " << (chains[i].t ? dfi(chains[i].t) : -1)<< "| ";
	   out << "type ";
	   switch(chains[i].type) {
	   case one: out << "1"; break;
	   case two_a: out << "2a"; break;
	   case two_b: out << "2b"; break;
	   case three_a: out << "3a"; break;
	   case three_b: out << "3b"; break;
	   case unmarked: out << "X"; break;
	   default : assert(false);
	   }
	   out << "}}\"";
	   if (chains[i].is_marked)
		   out << " color = red";
	   out << "]" << std::endl;
	   if (chains[i].parent!=NULL)
		   out << "node" << i << " -- " << "node" << chains[i].parent->number;
	   out << std::endl;
#else
	   out << "node" << i << " [label =\"{";
	   out << i << "|";
	   out << "{ s " << chains[i].s->id() << "} ";
	   out << "{ t " << chains[i].t->id() << "} ";
	   out << "{ type" << chains[i].type << "}";
	   out << "}\"]" <<std::endl;
	   out << "node" << i << " -- " << "node" << ((chains[i].parent == NULL) ? -1 : chains[i].parent->number);
	   out << std::endl;
#endif
   }

   for(unsigned int i = 0; i<num_chains; i++) {
	   caterpillar cp = caterpillars[i];
	   if (!cp.empty()) {
		   out << "subgraph clusterG" << i << " { " << std::endl;
		   for(caterpillar::item ch = cp.first(); ch != NULL; ch = cp.succ(ch)) {
			   out << "node" << cp.contents(ch)->number << std::endl;
		   }
		   out << "label = \"" << i << "\"}" << std::endl;
	   }

   }

   out << "}" << std::endl;
}


bool schmidt_is_triconnected(ugraph& g) {
	schmidt_triconnectivity t(g);
	std::fstream f("./blablubb.dot", std::ios::out);
	t.dfs_tree_to_dot(f);
	std::fstream f2("./chains.dot", std::ios::out);
	t.chain_tree_to_dot(f2);
	return true;
}
