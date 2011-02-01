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
#include <vector>

using namespace leda;

class schmidt_triconnectivity {


	edge_array<int> belongs_to_chain;
	edge_array<bool> is_frond;

	node_array<unsigned int> dfis;
	node* const node_at; //dfi reverse map
	node_array<edge> parent;
	node_array<int> inner_of_chain; // only inner nodes of chains are marked.
	node_array<slist<chain*>* > type_3;
	node_array<bool> is_real;

	std::vector<chain*> chains;
	caterpillar* const caterpillars;
	slist<const chain* const > construction_sequence;
	ugraph& the_graph;

public:
	schmidt_triconnectivity(ugraph& graph);

	~schmidt_triconnectivity(void);

	certificate* certify(void);

	void partition_into_segments(const chain* current_chain, const slist<chain*>& type3, h_array<unsigned int, slist<chain*> >& segments, h_array<unsigned int, slist<node> >& attachment_vertices, const int_set& set_children12, h_array<unsigned int,bool>& intersection_12_free);

	void add_easy_segments(const chain* current_chain, slist<chain*>& type3, const h_array<unsigned int, bool>& intersection_12_free);

	void add_hard_segments(const chain* current_chain, h_array<unsigned int, slist<node> >& attachment_vertices, slist<chain*>& children12);

	unsigned int dfi(const node v) const;
	void set_dfi(const node v, const unsigned int d);
	node parent_node(const node v) const;
	bool in_subdivision(const chain& c) const;
	bool in_subdivision(const chain* c) const;
	bool in_subdivision(const node n) const;

	void chain_decomposition(void);

	/* In this function we do a DFS on the graph to
	 * * calculate DFIs for each vertex
	 * * finding an initial K32 subdivision
	 */
	void initial_dfs(void);

	/* takes the first inner vertex of a chain (except C0?) and walks upwards in the tree until the edge to the parent is contained in another chain
	 * Marks edges on the path as belonging to the current chain. Marks inner nodes as inner nodes of the chain. Sets chain.t */
	void mark_path(const node start, const unsigned int the_chain);

	inline bool contained_in_chain(const node v, const chain* chain);

	chain_type classify_chain(const unsigned int chain_number);

	bool find_a_cycle_to_root(node& start_node, edge& backedge, unsigned int& number_seen, edge_array<bool>& seen_edge, bool continue_after_found, const unsigned int current_chain);

	void add_to_subdivision(chain* c);

   void dfs_tree_to_dot(std::ostream& out);

   void chain_tree_to_dot(std::ostream& out);

   void create_chain(edge e, unsigned int chain_number);

   void mark_as_frond(const edge e);
};

#endif /* SCHMIDT_HPP_ */
