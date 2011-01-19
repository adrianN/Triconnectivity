/*
 * InterlockingIntervals.hpp
 *
 *  Created on: Jan 14, 2011
 *      Author: adrian
 */

#ifndef INTERLOCKINGINTERVALS_HPP_
#define INTERLOCKINGINTERVALS_HPP_

#include "LEDA/graph/ugraph.h"
#include "LEDA/core/stack.h"
#include "LEDA/core/list.h"
#include "dfs.hpp"
#include "utilities.hpp"
#include <vector>
using namespace std;
using namespace leda;


template <typename A> class interval{
public:
	A cont;
	unsigned int bounds[2];
	leda::node represented_by[2];
	interval(unsigned int b[], A content) : cont(content) {
		std::cout << " create interval " << b[0] << " " << b[1] << " " << content << std::endl;
		assert(b[0]<b[1]);
		bounds[0] = b[0];
		bounds[1] = b[1];
		represented_by[0] = NULL;
		represented_by[1] = NULL;
	}
};

template<typename A> class Order {
public:

	static bool compute_order(const slist<interval<A>*>& sorted_asc,
							  const slist<interval<A>*>& sorted_dsc,
							  const slist<slist<interval<A>*>* >& equivalent_intervals,
							  interval<A>* start,
							  std::vector<interval<A>* >& output)
	{
		ugraph g;
		connect_forest(g, sorted_asc, 1);

		connect_forest(g, sorted_dsc, 0);

		glue_trees(g, sorted_asc);

		glue_equivalence_classes(g, equivalent_intervals);

		node* ordered = new node[g.number_of_nodes()];
		const node startnode = start->represented_by[0];
		bool is_connected = dfs_order(g, startnode, ordered);
		if (!is_connected) {
			delete[] ordered;
			return false;
		}

		node_array<interval<A>*> belongs_to(g, NULL);
		{
			interval<A>* i;
			forall(i, sorted_asc) {
				belongs_to[i->represented_by[0]] = i;
			}
		}


		for(unsigned int i=0; i<(unsigned int)g.number_of_nodes(); i++) {
			output.push_back(belongs_to[ordered[i]]);
		}
		delete[] ordered;
		return true;
	}

private:
	static node_array<A> labels(const ugraph& g, const slist<interval<A>*> list, unsigned int i) {
		node_array<A> l(g);
		interval<A>* val;
		forall(val, list) {
			assert(val->represented_by[i]!=NULL);
			l[val->represented_by[i]] = val->cont;
		}
		return l;
	}

	//returns an array with node->interval assocs
	static void connect_forest(ugraph& g, const slist<interval<A>* >& input_list, unsigned int side) {
		std::vector<interval<A>* > input_array;
		{	// create a node for each interval and copy the intervals into a vector
			interval<A>* val=NULL;
			forall(val, input_list) {
				assert(val!=NULL);
				input_array.push_back(val);
				node node_for_interval = g.new_node();
				assert(val->represented_by[side] == NULL);
				val->represented_by[side] = node_for_interval;
			}
		}

		//solve ANLV problem and build forest

		stack<interval<A>* > s;
		for(int i=input_array.size()-1; i>=0; i--) {
			interval<A>* cur_interval = input_array[i];
			assert(cur_interval!=NULL);
			while(!s.empty() &&	((side == 1 && s.top()->bounds[side] <= cur_interval->bounds[side]) || (side == 0 && s.top()->bounds[side] >= cur_interval->bounds[side])) ) { // < ?
				s.pop();
			}

			if (!s.empty() && ((side == 1 && s.top()->bounds[(side+1)%2] < cur_interval->bounds[side]) || (side == 0 && s.top()->bounds[(side+1)%2] > cur_interval->bounds[side]))) {
				g.new_edge(cur_interval->represented_by[side], s.top()->represented_by[side]); //top becomes parent
			}

			s.push(cur_interval);
		}
	}

	static void glue_trees(ugraph& g, const slist<interval<A> *> input) {
		interval<A>* interval;
		forall(interval, input) {
			interval->represented_by[0] = merge_nodes(g,interval->represented_by[0], interval->represented_by[1]);
			interval->represented_by[1] = NULL;
		}
	}

	static void glue_equivalence_classes(ugraph& g, const slist<slist<interval<A>*>* > equivalent_intervals) {
		slist<interval<A>* >* a_class;
		forall(a_class,equivalent_intervals) {
			slist<node> nodes_to_merge;
			interval<A>* i;
			forall(i, *a_class) {
				assert(i->represented_by[1]==NULL);
				nodes_to_merge.append(i->represented_by[0]);
			}
			node merged_node = merge_nodes(g, nodes_to_merge);
			forall(i,*a_class) {
				i->represented_by[0] = merged_node;
			}
		}

	}
};

#endif /* INTERLOCKINGINTERVALS_HPP_ */
