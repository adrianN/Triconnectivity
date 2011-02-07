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
	is_real(graph,false),
	caterpillars(new caterpillar[graph.number_of_edges() - graph.number_of_nodes() +2]), //TODO not space efficient, at most #fronds/2 caterpillars
	the_graph(graph),
	cert(new certificate(graph,this))
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


	assert(is_connected(the_graph));
	assert(the_graph.number_of_nodes()>=4);
	initial_dfs();
	chain_decomposition();

	//Assert property A TODO handle without asserts
	const unsigned int number_chains = chains.size();
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
		for(unsigned int i = 0; i<number_chains; i++) {
			if (chains[i]->get_s() == chains[i]->get_t()) {
				std::cout << "chain " << i << " forms a cycle" << std::endl;
				num_cycles++;
			}
		}
		assert(num_cycles==0 && "not biconnected"); //should already happen during the decomposition. Not 1 as in thesis, cause C0
		break;
	}
	default: break;
	}

	certify();
}

schmidt_triconnectivity::~schmidt_triconnectivity(void) {
	for(unsigned int i=0; i<chains.size(); i++)
		delete chains[i];
	delete[] caterpillars;
	delete[] node_at;
}

auto_ptr<certificate> schmidt_triconnectivity::certify(void) {
	std::cout << "\nStarting proper certification process" << std::endl;
	for(unsigned int i=0; i<(unsigned int)chains.size(); i++) {
		chain* current_chain = chains[i];
		std::cout << "Considering the children of chain " << current_chain->number << " (" << current_chain->children12.size() << ")"<< std::endl;
		assert(in_subdivision(current_chain));

		//Compute children 12, the set of all children of type 1 or 2 of current_chain that are not yet added to the subdivision
		slist<chain*> children12;
		int_set set_children12(chains.size());
		{	chain* child;
			forall(child,current_chain->children12) {

				if (!in_subdivision(child)) {
					switch (child->type) {
					case two_a:
						/* preemptively add chains of type 2a with real t(C). As they are backedges with parent in S they have prop. 30a. As t(C)
						 * is real and a proper descendant of t(current_chain) (which is real) by definition of type 2 chains, they have property 30b.
						 */
						if (is_real[child->get_t()]) {
							add_with_ancestors(child);
							continue;
						}
					case one: case two_b:
						children12.append(child);
						set_children12.insert(child->number); // TODO also for 2a that are added to subdivision?
						break;
					default: break;
					}
				}
			}
		}

		//compute type3, the set of chains of type 3 that start at an inner vertex of current_chain
		slist<chain*>& type3 = current_chain->type3;
		chain_node_iterator it(current_chain,this);
//		for(node v = it.next(); v!=NULL ; v = it.next()) { //<=?
//			if (type_3[v] != NULL)
//				type3.conc(*type_3[v]); //TODO order?
//		}
//		//this can probably be avoided by better traversal
//		type3 = bucket_sort<chain*,asc>(type3,get_number,0,chains.size());//this is bad, takes linear time! TODO

		{ chain* c; std::cout << "type 3 chains " ;
			forall(c,type3) {
				std::cout << c->number << " ";
			}
			std::cout << std::endl;
			std::cout << "Children12 ";
			forall(c,children12) {
				std::cout << c->number << " ";
			}
			std::cout << std::endl;
		}


		// Tag each chain in type3 with the segment it's contained in.
		// Compute the attachment vertices of each hard segment
		// The attachment vertices are the end vertices of the minimal chain and of the chains in type3 that are in the segment
		h_array<unsigned int, slist<node> > attachment_vertices;
		h_array<unsigned int, slist<chain*> > segment_chains;
		h_array<chain*, unsigned int> segment;

		partition_into_segments(current_chain, type3, attachment_vertices, segment_chains, segment, set_children12);

		//Add the easy segments, i.e. those without chains from children12
		add_easy_segments(current_chain, type3, segment, set_children12);

		// Add the hard clusters
		// we take the attachment vertices and generate a set of intervals for each segment. Then we compute an order
		// for the intervals such that we can add the minimal chain (and thus all others in the segment) such that there is at least one real vertex between
		// the end nodes of the chain.

		add_hard_segments(current_chain, attachment_vertices, segment_chains, children12);


	}
	return cert;
}

void schmidt_triconnectivity::partition_into_segments(const chain* current_chain, const slist<chain*>& type3, h_array<unsigned int, slist<node> >& attachment_vertices, h_array<unsigned int, slist<chain*> >& segment_chains, h_array<chain*,unsigned int>& segment, const int_set& children12) {
	node_array<int> minimal_chain(the_graph,-1);
	{	node n; //leda bugfix
		forall_nodes(n,the_graph) {
			minimal_chain[n] = -1;
		}
	}

	std::cout << "partitioning into segments" << std::endl;
	chain* t3_chain;
	forall(t3_chain, type3) {
		std::cout << "t3chain " << t3_chain->number << std::endl;
		//Find the minimal chain

		node v = t3_chain->get_t();
		while (!in_subdivision(parent_node(v)) && minimal_chain[v] < 0) {
			v = parent_node(v);
		};

		const bool new_segment = minimal_chain[v] < 0;
		if (new_segment) {
			minimal_chain[v] = inner_of_chain[v];
		}

		assert(minimal_chain[v] >= 0);
		const unsigned int m_chain = (unsigned int)minimal_chain[v];
		const bool hard = children12.member(m_chain);

		std::cout << "\tminchain " << m_chain << " (" << ((hard)? "hard" : "easy") << ")" << std::endl;

		//mark all nodes on the path to the minimal chain with the minimal chain to get linear running time.
		v = t3_chain->get_t();
		while (!in_subdivision(parent_node(v)) && minimal_chain[v] < 0) {
			minimal_chain[v] = m_chain;
			v = parent_node(v);
		}

		segment[t3_chain] = m_chain;
		segment_chains[m_chain].append(t3_chain);

		if (hard) { //only compute the attachment vertices of hard segments. Easy segments can be added without
			if (new_segment) { // new segment, add attachment vertices of minimal chain
				const chain* min_chain = chains[m_chain];
				assert(contained_in_chain(chains[m_chain]->get_s(), current_chain) && "m_chain can't be in children12");
				assert(contained_in_chain(chains[m_chain]->get_t(), current_chain) && "m_chain can't be in children12");
				attachment_vertices[m_chain].append(min_chain->get_s());
				attachment_vertices[m_chain].append(min_chain->get_t());
			}

			assert(contained_in_chain(chains[m_chain]->get_s(), current_chain) && "Start vertex of type3 chain is not an attachment vertex, it can't be in type3");
			attachment_vertices[m_chain].append(t3_chain->get_s());
		}
	}

}

void schmidt_triconnectivity::add_easy_segments(const chain* current_chain, slist<chain*>& type3, h_array<chain*,unsigned int> const & segment, const int_set& children12) {
	std::cout << "Adding easy segments" << std::endl;
	slist<chain*>::item prev_it=NULL;
	for(slist<chain*>::item cur_it = type3.first(); cur_it!=NULL; (prev_it=cur_it, cur_it=cur_it == NULL? type3.first() : type3.next_item(cur_it))){
		chain* t3_chain = type3.contents(cur_it);
		std::cout << "Considering chain " << t3_chain->number << std::endl;

		if (children12.member(segment[t3_chain])) continue; //not a hard segment

		add_with_ancestors(t3_chain);

		if (prev_it!=NULL) {
			std::cout <<  "Deleting chain\n" << std::endl;
			type3.del_succ_item(prev_it);
			cur_it=prev_it;
		} else {
			std::cout << "Popping chain\n" << std::endl;
			type3.pop();
			cur_it = NULL;
		}
	}
}

void schmidt_triconnectivity::add_with_ancestors(chain* a_chain) {
	decompose_to_bg_paths(a_chain);

	//add all ancestors of t3_chain that are still in H, in order of <
	stack<chain*> ancestors_in_segment;
	for(chain* cur_chain=a_chain; cur_chain->get_parent() != NULL && !cur_chain->in_subdivision ; cur_chain = cur_chain->get_parent()) {
		//assert(cur_chain->number > segment[t3_chain]);
		ancestors_in_segment.push(cur_chain);
	}

	while(ancestors_in_segment.size()>0) {
		add_to_subdivision(ancestors_in_segment.pop());
	}
}

void schmidt_triconnectivity::decompose_to_bg_paths(const chain* a_chain) {
	switch(a_chain->type) {
	case two_b: /* these must be part of caterpillars and get decomposed */ break;
	case three_b: {
		caterpillar& p = caterpillars[a_chain->number];
		assert(p.parent!=NULL);

		std::cout << "caterpillar with parent " << p.parent << std::endl;
		//test which kind of good caterpillar Ci we have
		enum caterpillartype {
			type1, 	//type 1: s(Ci) is on t(Ck) ->T s(Ck)
					//same as dfi(s(ck)) < dfi(sCi) < dfi(tck)
			type2	//type 2: there is a real vertex on s(Ci) ->ck s(Ck)
		} type;

		if (!p.parent->in_subdivision) {
			std::cout << "caterpillar is not good" << std::endl;
			assert(false && "graph is not triconnected");
		} else 	if (dfi(p.parent->get_s()) < dfi(a_chain->get_s()) && dfi(a_chain->get_s()) < dfi(p.parent->get_t())) {
			std::cout << "this is a type 1 caterpillar" << std::endl;
			type = type1;
		} else {
			chain_node_iterator it(p.parent, this);
			node n;
			it.next(); // discard s(ck)
			for(n = it.next(); !(n == NULL || is_real[n] || n == a_chain->get_s()); n= it.next()) { /* nothing */	}
			std::cout << "n: " << dfi(n) << " " << is_real[n] << " " << dfi(a_chain->get_s())<< std::endl;
			if (n == NULL) {
				std::cout << __LINE__ <<  " a bad caterpillar" << std::endl;
				assert(false && "graph is not triconnected");
			}
			else if (is_real[n]) {
				type = type2;
				std::cout << "good caterpillar of type 2" << std::endl;
			} else {
				std::cout << __LINE__ <<  " a bad caterpillar" << std::endl;
				assert(false && "graph is not triconnected");
			}
		}

		// fp: Ci + tree path up to first node in parent of caterpillar
		list<edge> ci; //edges in the chain ci
		list<edge> tci_y; //tree path from t to y
		h_array<edge,bool> on_tci_y(false);
		for(edge e = a_chain->first_edge; !contained_in_chain(source(e),p.parent) || dfi(target(e)) >= dfi(a_chain->get_t()); e = parent[target(e)]) {
			if (dfi(target(e))>=dfi(a_chain->get_t())) { // we're still on ci
				ci.append(e);
			} else {
				tci_y.append(e);
				on_tci_y[e] = true;

			}

		}

		list<edge> d0; //the parent of ci (not the caterpillar), up to tci
		{
			const chain* parent_chain = a_chain->get_parent();
			for(edge e = parent_chain->first_edge; source(e) != a_chain->get_t(); e = parent[target(e)]) {
				std::cout << "pc " << dfi(source(e)) << "  " << dfi(target(e)) << std::endl;
				d0.append(e);
			}
		}

		switch(type) {
		case type1:{
			std::cout << "type 1 caterpillar" << std::endl;
			ci.conc(tci_y, leda::behind);
			cert->add_bg_path(ci);
			cert->add_bg_path(d0);

		} break;
		case type2:{
			std::cout << "type 2 caterpillar" << std::endl;

			//d0.reverse();
			ci.reverse();
			d0.conc(ci, leda::behind);
			cert->add_bg_path(d0);
			cert->add_bg_path(tci_y);
		}
		}

		//now it is ok to add the parent chains. They are stored in the caterpillar.
		p.pop(); //discard the immediate parent of a_chain, it was processed as d0
		chain* c;
		forall(c,p) {
			//we only want the part of the chain until it hits the path from tci->y
			list<edge> di;
			for(edge e = c->first_edge; !on_tci_y[e]; e = parent[target(e)]) {
				di.append(e);
			}
			cert->add_bg_path(di);
		}


	} break;

	case one: case two_a: case three_a: case unmarked: cert->add_bg_path(a_chain);
	}
}


void schmidt_triconnectivity::add_hard_segments(const chain* current_chain, h_array<unsigned int, slist<node> > const & attachment_vertices, h_array<unsigned int, slist<chain*> > const & segment_chains, slist<chain*>& children12) {
	//we start by computing an order of the segments
	// we map the nodes of the current_chain to ints from 1 to |current_chain|
	// we add an artificial interval for each real node on the chain
	std::cout << "Adding hard segments" << std::endl;
	assert(current_chain->in_subdivision);
	if (children12.size() == 0 && attachment_vertices.size() == 0) {
		std::cout << "No hard segments" << std::endl;
		return;
	}

	node_array<unsigned int> mapping(the_graph,0);
	//std::vector<node> reverse_mapping;
	//reverse_mapping.push_back(NULL); //mapping starts at 1
	unsigned int vertices_in_chain = 0;
	{
		chain_node_iterator it(current_chain, this);
		//TODO change this to properly compute the dependent path
		std::cout << "Mapping and artificial intervals\n";
		for(node n=it.next(); n!=NULL; n=it.next()) {
			vertices_in_chain++;
			mapping[n] = vertices_in_chain;
			std::cout << "map " << dfi(n) << " -> " << vertices_in_chain << std::endl;
			//reverse_mapping.push_back(n);
		}
	}
	//the endpoints of the dependent path
	unsigned int min_n = vertices_in_chain+1, max_n = 0;
	slist<interval<chain* > * > intervals;
	std::cout << "Intervals for remaining children12" << std::endl;
	{ chain* c;
		forall(c,children12) {
			assert(!c->in_subdivision);
			if (attachment_vertices.defined(c->number)) continue; // avoid double-adding chains from children12 that are minchains of hard segments
			std::cout << c->number << std::endl;
			unsigned int m1,m2;
			m1 = mapping[c->get_s()];
			m2 = mapping[c->get_t()];
			min_n = std::min(min_n,std::min(m1,m2));
			max_n = std::max(max_n,std::max(m1,m2));
			intervals.append(new interval<chain*>(m1, m2,c));
		}
	}

	//iterate over the attachment vertices of all chains and generate the intervals
	//we generate one interval from the lowest attachment vertex to all others and one from the highest to all others
	std::cout << "Intervals for the attachment vertices of type 3 segments \nThere are " << attachment_vertices.size() << " such segments" << std::endl;
	slist<slist<interval<chain*>* >* > equivalent_intervals;
	{	unsigned int m_chain;
		forall_defined(m_chain, attachment_vertices) {

			if (chains[m_chain]->in_subdivision) // this was an easy chain and has already been added.
				continue;
			std::cout << "Attachment vertices for " << m_chain << std::endl;
			slist<node> const & vertices = attachment_vertices[m_chain];
			slist<unsigned int> mapped_vertices;
			node n;
			//TODO this can be simplified by walking in the right direction and avoiding one bucket sort and the intermediate list
			forall(n,vertices) {
				std::cout << " " << mapping[n];
				mapped_vertices.append(mapping[n]);
			}

			std::cout << std::endl;

			{
				unsigned int (*id)(const unsigned int& a) = (identity<unsigned int>);
				mapped_vertices = bucket_sort<unsigned int,asc>(mapped_vertices, id, 0, vertices_in_chain);
			}

			slist<interval<chain*>*>* equivalence_class = new slist<interval<chain*>*>();
			equivalent_intervals.append(equivalence_class);

			//generate the intervals
			const unsigned int smallest = mapped_vertices.contents(mapped_vertices.first());
			const unsigned int largest = mapped_vertices.contents(mapped_vertices.last());

			min_n = std::min(min_n,smallest);
			max_n = std::max(max_n,largest);

			interval<chain*>* val = new interval<chain*>(smallest,largest,chains[m_chain]);
			intervals.append(val);
			equivalence_class->append(val);

			slist<unsigned int>::item it = mapped_vertices.first();
			it = mapped_vertices.succ(it);

			for(unsigned int i=1; i <= (unsigned int)mapped_vertices.size()-2; i++) {
				val = new interval<chain*>(smallest, mapped_vertices.contents(it), chains[m_chain]);
				intervals.append(val);
				equivalence_class->append(val);

				val = new interval<chain*>(mapped_vertices.contents(it), largest, chains[m_chain]);
				intervals.append(val);
				equivalence_class->append(val);
			}
		}
	}


	{	interval<chain*>* start_interval = NULL;
		//create artificial intervals for vertices on the dependend path that are already real
		std::cout << "real vertices on the dependent path " << min_n << " " << max_n << std::endl;
		slist<interval<chain*>*>* equivalence_class = new slist<interval<chain*>* >();
		{ 	chain_node_iterator it(current_chain,this);
			for(node n = it.next(); n!=NULL; n=it.next()) {
				if (is_real[n] && mapping[n] >= min_n && mapping[n] <=max_n) {
					interval<chain*>* val  = new interval<chain*>(0,mapping[n], NULL); // artificial intervals for real vertices start at 0
					intervals.append(val);
					equivalence_class->append(val);
				}
			}
			start_interval = equivalence_class->contents(equivalence_class->last());
			// the artificial vertices should also be in an equivalence class, as it is ok if some of them don't overlap other intervals.
			equivalent_intervals.append(equivalence_class);
		}

		// sort the intervals
		{
			unsigned int (*fc)(interval<chain*>* const &) = &interval<chain*>::first_component;
			unsigned int (*sc)(interval<chain*>* const &) = &interval<chain*>::second_component;

			/* because these intervals don't have distinct endpoints as assumed in the paper, we need a more complicated ordering. Luckily this
			 * can also be achieved using bucket sort. We just do two passes over the intervals, such that if we for example sort primarily by fc
			 * asc, we also order by sc dsc. This way if L1 = L2 but R1 < R2, R2 comes first in the ordering. As described in "Fast algorithm for interlocking sets (1989)"f
			 */
			slist<interval<chain*>*> intervals_asc = bucket_sort<interval<chain*>*,asc, dsc>(intervals,fc,sc,0,vertices_in_chain);
			slist<interval<chain*>*> intervals_dsc = bucket_sort<interval<chain*>*,dsc, asc>(intervals,sc,fc,0,vertices_in_chain);

		// order the intervals in fancy order and add them
			std::vector<interval<chain*>*> ordered;
			bool order_exists = Order<chain*>::compute_order(intervals_asc, intervals_dsc, equivalent_intervals, start_interval, ordered);
			assert(order_exists); //other graph not triconnected.
			std::cout << "order is: ";
			for(unsigned int i=0; i<ordered.size(); i++) {
				std::cout << ordered[i]->cont << " ";
			}
			std::cout << std::endl;
			for(unsigned int i = 0; i< ordered.size(); i++) {
				if (ordered[i]->cont!=NULL) { // not an artificial interval
					std::cout << ordered[i] << std::endl;

					add_with_ancestors(ordered[i]->cont);
					if (segment_chains.defined(ordered[i]->cont->number)) {
						chain* chain_in_segment;
						forall(chain_in_segment,segment_chains[ordered[i]->cont->number]) {
							add_with_ancestors(chain_in_segment);
						}
					}
					//add_to_subdivision(ordered[i]->cont);
				}
			}
		}
	}
	// free resources
	{
		interval<chain*>* val;
		forall(val, intervals) {
			delete val;
		}
		slist<interval<chain*>*>* equivalence_class;
		forall(equivalence_class, equivalent_intervals) {
			delete equivalence_class;
		}
	}


}

unsigned int schmidt_triconnectivity::dfi(const node v) const {assert(v!=NULL); return dfis[v];}
void schmidt_triconnectivity::set_dfi(const node v, const unsigned int d) { assert(v!=NULL); assert(dfis[v] == 0); dfis[v] = d; node_at[d] = v;}

node schmidt_triconnectivity::parent_node(const node v) const {
	assert(v!=NULL);
	edge e = parent[v];
	if (e!=NULL) {
		assert(!is_frond[e]);
		assert(source(e) == v);
		return target(e);
	} else {

		return NULL;
	}
}

bool schmidt_triconnectivity::in_subdivision(const chain& c) const {
	return c.in_subdivision;
}

bool schmidt_triconnectivity::in_subdivision(const chain* c) const {
	return c->in_subdivision;
}

bool schmidt_triconnectivity::in_subdivision(const node n) const {
	return dfi(n)==1 || in_subdivision(chains[inner_of_chain[n]]); //every node is the inner node of some chain, except maybe the root. I'm not sure about that TODO
}

void schmidt_triconnectivity::create_chain(edge e, unsigned int chain_number) {
	assert(is_frond[e]);
	chain* current_chain = new chain();
	chains.push_back(current_chain);
	assert(chains.size()==chain_number+1);

	//fronds go from up to down, hence chain.s = source, chain.t = target
	std::cout << "edge " << dfi(source(e)) << "->" << dfi(target(e)) << " belongs to " << chain_number << "\n";
	assert(belongs_to_chain[e]<0);
	belongs_to_chain[e] = chain_number;
	current_chain->set_s(source(e));

	std::cout << "create chain from frond " << dfi(source(e)) << "->" << dfi(target(e)) << std::endl;

	current_chain->first_edge = e;

	mark_path(target(e),chain_number); //sets .t

	assert(current_chain->get_s()!=current_chain->get_t() && "graph is not biconnected");
	assert(dfi(current_chain->get_s()) < dfi(current_chain->get_t()));

	switch (classify_chain(chain_number)) {
	case three_a: case three_b: {
		const node s = current_chain->get_s();
		chain* chain_at_s = chains[inner_of_chain[s]];
		chain_at_s->type3.append(current_chain);
	}

	default: break;
	}
}

void schmidt_triconnectivity::chain_decomposition(void) {
	//Calculate the list of fronds in DF order of their starting vertex and find the chain for each frond. This can probably be incorporated into the DFS.
	unsigned int chain_number = 3; //the first three chains are found during dfs

	for(unsigned int cur_dfi=1; cur_dfi<(unsigned int)the_graph.number_of_nodes()+1; cur_dfi++) {
		node current_node = node_at[cur_dfi];
		assert(current_node!=NULL);
		edge e;
		forall_adj_edges(e,current_node) {
			if (!is_frond[e] || belongs_to_chain[e] >= 0) // the first three chains already exist
				continue;

			create_chain(e,chain_number);

			assert(chains[chain_number]!=NULL);
			assert(chains[chain_number]->get_s()!=NULL);
			assert(chains[chain_number]->get_t()!=NULL);
			assert(dfi(chains[chain_number]->get_s()) <= dfi(chains[chain_number]->get_t())); // == if not triconnected.
			assert(chains[chain_number]->number == chain_number);
			assert(cur_dfi==0 || chains[chain_number]->first_edge !=NULL);

			chain_number++;
		}
	}

	std::cout << "found " << chain_number-1 << " many chains" << endl;


	std::fstream f("./blubb.dot", std::ios::out);
	dfs_tree_to_dot(f);
	f.close();

	std::fstream g("./chains.dot", std::ios::out);
	chain_tree_to_dot(g);
	g.close();

#ifndef NDEBUG
	//Chain decompositition partitions edges
	{	edge e;
		forall_edges(e,the_graph) {
			assert(belongs_to_chain[e] >= 0 && (unsigned int)belongs_to_chain[e] < chains.size());
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
		inner_of_chain[root] = 0; //TODO is this so?
		assert(the_graph.degree(root)>=3);

		unsigned int number_seen=1;
		set_dfi(current_node, number_seen++);

		parent[current_node] = NULL;
		edge_array<bool> seen_edge(the_graph, false);



		/* first find the first fundamental cycle back to the root
		 * do a normal dfs
		 */
		std::cout << "first cycle " << std::endl;
		edge chain1_first_edge;
		if (!find_a_cycle_to_root(current_node, chain1_first_edge, number_seen, seen_edge, false,  1)) {
			assert(false && "graph is not triconnected");
		}

		assert(current_node == target(chain1_first_edge));

		std::cout << "node a: "  << current_node->id() << std::endl;
		/* When the first cycle is found we start a new DFS for every node on the path from a upwards to the root. When one
		 * of those searches encounters another backedge to the root we know that the starting node is the LCA.
		 *
		 * TODO we can also start labeling edges as belonging to a path
		 */

		node inner_search_cur_node;
		node lca = NULL;
		node node_b = NULL;
		edge chain2_first_edge = NULL;
		while(current_node != NULL) { //todo switch from node to edge
			std::cout << "walk up " << current_node->id() <<std::endl;


			inner_search_cur_node = current_node;
			edge e;
			if (find_a_cycle_to_root(inner_search_cur_node,e, number_seen, seen_edge, true, 2) && lca == NULL && node_b == NULL) {
				std::cout << "found a cycle " << std::endl;
				lca = current_node;
				node_b = inner_search_cur_node;
				chain2_first_edge = e;
				assert(target(chain2_first_edge)==node_b);

			}

			current_node = parent_node(current_node);

		}
		assert(lca!=NULL);
		assert(node_b!=NULL);
		assert(target(chain2_first_edge)==node_b);


		std::cout << "Root " << root->id() << std::endl;
		std::cout << "LCA " << lca->id() << std::endl;
		std::cout << "Chain one fe " << dfi(source(chain1_first_edge)) << " " << dfi(target(chain1_first_edge)) << std::endl;
		std::cout << "Chain two fe " << dfi(source(chain2_first_edge)) << " " << dfi(target(chain2_first_edge)) << std::endl;

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
		chains.push_back(new chain());

		chains[0]->number = 0;
		chains[0]->set_parent(NULL);
		chains[0]->first_edge = parent[lca];


		//TODO, do we need these at all? C0 stays unclassified
		chains[0]->set_s(lca);
		mark_path(lca, 0);

		add_to_subdivision(chains[0]);
		cert->add_bg_path(chains[0]);

		create_chain(chain1_first_edge,1);
		create_chain(chain2_first_edge,2);

		add_to_subdivision(chains[1]);
		add_to_subdivision(chains[2]);
		cert->add_bg_path(chains[1]);
		cert->add_bg_path(chains[2]);

		assert(chains[0]->get_s() == chains[2]->get_t());
		assert(chains[0]->get_t() == root);
		assert(chains[1]->get_s() == root);
		assert(chains[2]->get_s() == root);
		assert(chains[1]->get_t() == lca);
		assert(chains[2]->get_t() == lca);


		for(unsigned int i= 0 ; i<3; i++)
			classify_chain(i);
}

/* takes the first inner vertex of a chain (except C0?) and walks upwards in the tree until the edge to the parent is contained in another chain
 * Marks edges on the path as belonging to the current chain. Marks inner nodes as inner nodes of the chain. Sets chain.t */
void schmidt_triconnectivity::mark_path(const node start, const unsigned int the_chain) {
	chains[the_chain]->number = the_chain;
	std::cout << "Chain " << the_chain << "\n\t";
	node current_node = start;
	for(edge current_edge = parent[start]; ; (current_node = parent_node(current_node), current_edge = parent[current_node])) {

		if (current_edge == NULL || belongs_to_chain[current_edge]>=0) {
			//found an edge that belongs to a different chain, hence this chain ends here.
			assert((the_chain == 0) || ((unsigned int)belongs_to_chain[current_edge] < the_chain && "the parent of a chain must be a chain of smaller chain number"));
			chains[the_chain]->set_parent(current_edge == NULL ? NULL : chains[belongs_to_chain[current_edge]]);
			chains[the_chain]->set_t(current_node);
			assert(inner_of_chain[chains[the_chain]->get_t()]>=0);
			assert(the_chain == 0 || (unsigned int)inner_of_chain[chains[the_chain]->get_t()] == chains[the_chain]->get_parent()->number);
			break;
		}

		chains[the_chain]->is_backedge = false; // more than one edge contained

		assert(inner_of_chain[current_node] < 0 && "No node is inner node of >1 chains");
		inner_of_chain[current_node] = the_chain;
		//must be <0 due to if
		//std::cout << "edge " << dfi(source(current_edge)) << "->" << dfi(target(current_edge)) << " belongs to " << the_chain << "\n";

		assert(belongs_to_chain[current_edge]<0 && "No edge belongs to more than one chain");
		belongs_to_chain[current_edge] = the_chain;

	}

	std::cout << std::endl;
}

inline bool schmidt_triconnectivity::contained_in_chain(const node v, const chain* chain) {
	return (unsigned int)inner_of_chain[v] == chain->number || v == chain->get_s() || v == chain->get_t();
}

chain_type schmidt_triconnectivity::classify_chain(const unsigned int chain_number) {
	if (chain_number == 0)
		return unmarked; // first chain stays unmarked

	chain* const the_chain = chains[chain_number];
	chain* const the_parent = the_chain->get_parent();

	assert(contained_in_chain(the_chain->get_t(),the_parent));

	// t(the_chain) ->T s(the_chain) is contained in parent
	if (the_parent->number == 0 || (the_chain->get_s() != the_parent->get_s() && contained_in_chain(the_chain->get_s(), the_parent))) {
#ifndef NDEBUG
		//walk the path the_chain.t ->T the_chain.s and check if each node is inner_node in parent
		{
			node current_node = the_chain->get_t();
			assert((unsigned int)inner_of_chain[current_node] == the_parent->number);
			do {
				current_node = parent_node(current_node);
				assert(contained_in_chain(current_node,the_parent));
			} while (current_node != the_chain->get_s());
		}
#endif
		the_chain->type = one;
		return one;
	}

	//the chain starts at the same node as its parent, it is of type 2
	if (the_chain->get_s() == the_parent->get_s()) {
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
		//caterpillars[chain_number].append(the_chain);
		chain* c_jay = the_parent;
		while(c_jay->is_marked) {
			assert(c_jay->type = two_b);

			c_jay->is_marked = false;
			caterpillars[chain_number].append(c_jay);
			c_jay = c_jay->get_parent();
		}
		return three_b;
	}
}

/* do a dfs until an edge back to the root is encountered. if continue_after_found is true, continue anyway. */
bool schmidt_triconnectivity::find_a_cycle_to_root(node& start_node, edge& backedge, unsigned int& number_seen, edge_array<bool>& seen_edge, bool continue_after_found, const unsigned int current_chain) {
		stack<edge> first_cycle;

		/* push all neighbours of the current node */ //TODO is this necessary? Should be enough to push this node
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

				mark_as_frond(current_edge);

				if (dfi(next_node) == 1 && !found_cycle) { //backedge to the root
					start_node = current_node; //pass current node outside
					backedge = current_edge;
					assert(target(backedge) == start_node);
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
	assert(c!=NULL);
	std::cout << "Adding to subdivision " << c->number << std::endl;
	assert(!c->in_subdivision);
	assert((c->get_parent() == NULL && c->number==0) || c->get_parent()->in_subdivision); //modularity
	c->in_subdivision = true;
	is_real[c->get_s()] = true;
	is_real[c->get_t()] = true;

#ifndef NDEBUG //check that s and t really are real, i.e. have >=3 edges in the subdivision
	if (c->number > 2) {
	edge e;
	unsigned int count=0;
		forall_adj_edges(e,c->get_s()) {
			assert(belongs_to_chain[e]>=0 && (unsigned int)belongs_to_chain[e] < chains.size());
			if (chains[belongs_to_chain[e]]->in_subdivision) {
				count++;
			}
		}
		assert(count>=3 && "node doesn't have enough edges to become real");

		count=0;
		forall_adj_edges(e,c->get_t()) {
			assert(belongs_to_chain[e]>=0 && (unsigned int)belongs_to_chain[e] < chains.size());
			if (chains[belongs_to_chain[e]]->in_subdivision) {
				count++;
			}
		}
		assert(count>=3 && "node doesn't have enough edges to become real");
	}
#endif
}

void schmidt_triconnectivity::mark_as_frond(const edge e) {
	std::cout<< "new frond " << dfi(source(e)) << " "<< dfi(target(e)) << std::endl;
	assert(dfi(source(e))!=0 && "fronds have both ends explored");
	assert(dfi(target(e))!=0 && "fronds have both ends explored");
	assert(dfi(source(e))<dfi(target(e)) && "fronds go away from the root");
	is_frond[e] = true;
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
   const unsigned int num_chains = chains.size();
   for(unsigned int i = 0; i<num_chains; i++) {
#ifdef RECORD_NODES
	   out << "node" << i << " [shape = record, label =\"{";
	   out << i << "| {";
	   out << "s " << (chains[i]->get_s() ? dfi(chains[i]->get_s()) : -1) << "| ";
	   out << "t " << (chains[i]->get_t() ? dfi(chains[i]->get_t()) : -1)<< "| ";
	   out << "type " << chains[i]->type;
	   out << "}}\"";
	   if (chains[i]->is_marked)
		   out << " color = red";
	   out << "]" << std::endl;
	   if (chains[i]->get_parent()!=NULL)
		   out << "node" << i << " -- " << "node" << chains[i]->get_parent()->number;
	   out << std::endl;
#else
	   out << "node" << i << " [label =\"{";
	   out << i << "|";
	   out << "{ s " << chains[i]->get_s()->id() << "} ";
	   out << "{ t " << chains[i]->get_t()->id() << "} ";
	   out << "{ type" << chains[i]->type << "}";
	   out << "}\"]" <<std::endl;
	   out << "node" << i << " -- " << "node" << ((chains[i]->parent == NULL) ? -1 : chains[i].parent->number);
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

	return true;
}
