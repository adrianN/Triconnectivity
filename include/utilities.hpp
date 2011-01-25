/*
 * utilities.hpp
 *
 *  Created on: Oct 25, 2010
 *      Author: adrian
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include "LEDA/graph/ugraph.h"
#include "LEDA/core/slist.h"
#include <ostream>
#include <istream>
#include <memory>

void to_dot(const leda::ugraph& g, std::ostream& output);
void simple_to_dot(const leda::ugraph& g, std::ostream& output);

std::istream& operator>>(std::istream& input, leda::ugraph& graph);
leda::node merge_nodes(leda::ugraph& g, leda::node one, leda::node two);
leda::node merge_nodes(leda::ugraph& g, leda::slist<leda::node> nodes);
void glue_graphs(leda::ugraph& one, leda::ugraph& two, leda::node& merge_one, leda::node& merge_two);
void write_planar_code(leda::ugraph& g, std::ostream& out);
std::auto_ptr< leda::ugraph > triconnected_graph(const unsigned int n);
leda::slist<unsigned int> bucket_sort(leda::slist<unsigned int> elems, unsigned int start,  unsigned int end);
template<typename A> A identity(const A& a) { return a; }

template <typename A> leda::slist<A> bucket_sort(leda::slist<A> elems, unsigned int (*to_int)(const A&), unsigned int start,  unsigned int end) {
	A* buckets = new A[end-start+1];
	bool* filled = new bool[end-start+1];

	for(unsigned int i=0; i<(end-start+1); i++) {
		filled[i] = false;
	}
	A elem;
	forall(elem, elems) {
		const unsigned int elem_v = to_int(elem);
		assert(!filled[elem_v]);
		buckets[elem_v-start] = elem;
	}

	leda::slist<A> ret;
	for(unsigned int i=0; i<(end-start+1); i++) {
		if (filled[i])
			ret.append(buckets[start+i]);
	}
	delete[] buckets;
	delete[] filled;

	return ret;
}


template <typename A> void to_dot(const leda::ugraph& g, const leda::node_array<A>& labels, std::ostream& out) {

    out << "digraph G {" << std::endl;
    leda::node n;
    forall_nodes(n,g) {

        out << "\tnode" << n->id() << " [label = \"" << labels[n] << "\"]" << std::endl;
    }
    leda::edge e;
    forall_edges(e,g) {
        leda::node u = source(e);
        leda::node n = target(e);
        out << "\tnode" << n->id() << " -> " << "node" << u->id() << std::endl;
    }

    out << "}" << std::endl;
}
#endif /* UTILITIES_HPP_ */
