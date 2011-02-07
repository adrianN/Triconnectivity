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
bool graphs_isomorphic(leda::ugraph  & g1, leda::ugraph  & g2, leda::node_array<leda::node> const & map_2_to_1);
leda::edge smoothe(leda::ugraph& g, leda::node n); //returns the new edge
template<typename A> A identity(const A& a) { return a; }

enum ord { asc, dsc};

template <typename A, ord o> leda::slist<A> bucket_sort(leda::slist<A> elems, unsigned int (*to_int)(const A&), unsigned int start,  unsigned int end) {
	leda::slist<A>* buckets = new leda::slist<A>[end-start+1];
	A elem;
	forall(elem, elems) {
		const unsigned int elem_v = to_int(elem);
		buckets[elem_v-start].append(elem);
	}

	leda::slist<A> ret;
	switch(o) {
	case asc :
		for(unsigned int i=0; i<(end-start+1); i++) {
			ret.conc(buckets[start+i]);
		}
		break;
	case dsc:
		for(int i=end-start; i>=0; i--) {
			ret.conc(buckets[start+i]);
		}
		break;
	}
	delete[] buckets;
	assert(ret.size() == elems.size());
	return ret;
}


template <typename A, ord o1, ord o2> leda::slist<A> bucket_sort(leda::slist<A> elems, unsigned int (*to_int_first)(const A&), unsigned int (*to_int_second)(const A&), unsigned int start,  unsigned int end) {
	elems = bucket_sort<A,o2>(elems,to_int_second,start,end);
	return bucket_sort<A,o1>(elems,to_int_first,start,end);
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
