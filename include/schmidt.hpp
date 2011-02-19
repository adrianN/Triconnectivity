/*
 * schmidt.hpp
 *
 *  Created on: Jan 25, 2011
 *      Author: adrian
 */

#ifndef SCHMIDT_HPP_
#define SCHMIDT_HPP_

#include "LEDA/core/slist.h"
#include "LEDA/core/h_array.h"
#include "LEDA/core/int_set.h"
#include "LEDA/graph/ugraph.h"
#include "chain.hpp"
#include "caterpillar.hpp"
#include "certificate.hpp"
#include "not_triconnected_exception.hpp"
#include <vector>
#include <memory>

using namespace leda;
using namespace std;
class certificate; //NO IDEA why this is necessary.



class schmidt_triconnectivity {


	edge_array<int> belongs_to_chain;
	edge_array<bool> is_frond;

	node_array<unsigned int> dfis;
	node* node_at;
	node_array<edge> parent;
	node_array<int> inner_of_chain; // only inner nodes of chains are marked.
	node_array<bool> is_real;

	std::vector<chain*> chains;
	h_array<unsigned int, caterpillar> caterpillars;
	ugraph& the_graph;
	auto_ptr<certificate> cert;
	unsigned int chains_in_subdivision;

public:
	friend class chain_edge_iterator;
//	schmidt_triconnectivity(ugraph& graph) throw(not_triconnected_exception);
	schmidt_triconnectivity(ugraph& graph, node startnode) throw(not_triconnected_exception);

	~schmidt_triconnectivity(void);

	auto_ptr<certificate> certify(void);

	void partition_into_segments(
			const chain* current_chain,
			const list<chain*>& type3,
			h_array<unsigned int, slist<node> >& attachment_vertices,
			h_array<unsigned int, slist<chain*> >& segment_chains,
			h_array<chain*,unsigned int>& segment,
			const int_set& children12) throw();
	void add_with_ancestors(chain* a_chain) throw();
	void add_easy_segments(const list<chain*>& type3, h_array<chain*,unsigned int> const & segment, const int_set& children12);
	void decompose_to_bg_paths(const chain* a_chain) throw();

	void check_prop_b(void) const throw(not_triconnected_exception);

	void add_hard_segments(
			const chain* current_chain,
			h_array<unsigned int,
			slist<node> > const & attachment_vertices,
			h_array<unsigned int, slist<chain*> > const & segment_chains,
			slist<chain*>& children12) throw(not_triconnected_exception);

	unsigned int dfi(const node v) const throw();
	void set_dfi(const node v, const unsigned int d) throw();
	node parent_node(const node v) const throw();

	bool in_subdivision(node v) const throw();

	void chain_decomposition(void) throw(not_triconnected_exception);

	/* In this function we do a DFS on the graph to
	 * * calculate DFIs for each vertex
	 * * finding an initial K32 subdivision
	 */
	void initial_dfs(node startnode) throw(not_triconnected_exception);

	/* takes the first inner vertex of a chain (except C0?) and walks upwards in the tree until the edge to the parent is contained in another chain
	 * Marks edges on the path as belonging to the current chain. Marks inner nodes as inner nodes of the chain. Sets chain.t */
	void mark_path(const node start, const unsigned int the_chain) throw();

	bool contained_in_chain(const node v, const chain* chain) const throw();

	void classify_chain(const unsigned int chain_number) throw();

	bool find_a_cycle_to_root(node& start_node, edge& backedge, unsigned int& number_seen, edge_array<bool>& seen_edge, bool continue_after_found) throw();

	void add_to_subdivision(chain* c) throw();

   void dfs_tree_to_dot(std::ostream& out);

   void chain_tree_to_dot(std::ostream& out);

   void create_chain(edge e, unsigned int chain_number) throw(not_triconnected_exception);

   void mark_as_frond(const edge e) throw();
};

#endif /* SCHMIDT_HPP_ */
