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
#include "chain.hpp"
#include <vector>
#include <string>
#include <fstream>
using namespace std;
using namespace leda;




template <typename A> class interval{
public:
	A cont;
	unsigned int bounds[2];
	leda::node represented_by[2];
	interval(unsigned int l, unsigned int r, A content) : cont(content) {
		std::cout << "Fresh interval " << l << " " << r  << " " << cont << std::endl;
		if (l>r)
			std::swap(l,r);
		bounds[0] = l;
		bounds[1] = r;
		represented_by[0] = NULL;
		represented_by[1] = NULL;
	}
	interval(unsigned int b[], A content) : cont(content) {
		assert(b[0]<b[1]);
		bounds[0] = b[0];
		bounds[1] = b[1];
		represented_by[0] = NULL;
		represented_by[1] = NULL;
	}

	static unsigned int first_component(interval<A>* const & i) { return i->bounds[0]; }
	static unsigned int second_component(interval<A>* const & i) { return i->bounds[1]; }
};

template<typename A> std::ostream& operator<<(std::ostream& c, interval<A>* const& i) {
	if (i==NULL)
		c << "null";
	else
		c << "(" << i->bounds[0] << "," << i->bounds[1] << "," << i->cont << ")";
	return c;
}

template<typename A> class Order {

	public:

	static bool compute_order(const slist<interval<A>*>& sorted_asc,
							  const slist<interval<A>*>& sorted_dsc,
							  const slist<slist<interval<A>*>* >& equivalent_intervals,
							  const interval<A>* start,
							  std::vector<interval<A>* >& output)
	{static unsigned char b='0';
		assert(sorted_asc.size() == sorted_dsc.size());
		assert(sorted_asc.size() > 0);
		std::cout << "Computing overlap graph " << sorted_asc.size() << std::endl;
		ugraph g;
		connect_forest(g, sorted_asc, 1);
		{
			std::string s("a_forest");
			s += b;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}

		connect_forest(g, sorted_dsc, 0);
		{
			std::string s("b_forest");
			s += b;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}

		glue_trees(g, sorted_asc);
		{
			std::string s("a_graph");
			s += b;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}

		glue_equivalence_classes(g, equivalent_intervals);
		{
			std::string s("final_graph");
			s += b++;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}

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
	static node_array<interval<A>*> labels(const ugraph& g, const slist<interval<A>*> list) {
		node_array<interval<A>*> l(g);
		interval<A>* val;
		for(unsigned int i = 0; i<2; i++) {
			forall(val, list) {
				if(val->represented_by[i] != NULL)
					l[val->represented_by[i]] = val;
			}
		}
		return l;
	}

	static void connect_forest(ugraph& g, const slist<interval<A>* >& input_list, unsigned int side) {
		std::vector<interval<A>*> input_array;
		{	// create a node for each interval and copy the intervals into a vector
			interval<A>* val=NULL;
			//std::cout << "intervals: ";
			forall(val, input_list) {
				assert(val!=NULL);
				std::cout << val << " ";
				input_array.push_back(val);
				node node_for_interval = g.new_node();
				assert(val->represented_by[side] == NULL);
				val->represented_by[side] = node_for_interval;
			}
			std::cout << std::endl;
		}

		//solve ANLV problem and build forest

		stack<interval<A>* > s;
		for(int i=input_array.size()-1; i>=0; i--) {
			interval<A>* cur_interval = input_array[i];
			assert(cur_interval!=NULL);
			//fancy conditions to switch between <= and >= for different sides
			while(!s.empty() &&	((side == 1 && s.top()->bounds[1] <= cur_interval->bounds[1]) || (side == 0 && s.top()->bounds[0] >= cur_interval->bounds[0])) ) { // < ?
				s.pop();
			}

			//two intervals t, c overlap if (side 0) t.a < c.a < t.b < c.b or (side 1) c.a < t.a < c.b < t.b
			//by while loop top element has
			// 1: c.b < t.b
			assert(s.empty() || side != 1 || s.top()->bounds[1] > cur_interval->bounds[1]);
			// 0: t.a < c.a
			assert(s.empty() || side != 0 || s.top()->bounds[0] < cur_interval->bounds[1]);

			const bool s1_ta_le_cb = s.empty() || (side == 1 && s.top()->bounds[0] < cur_interval->bounds[1]); //t.a<c.b
			const bool s1_ca_le_ta = s.empty() || (side == 1 && cur_interval->bounds[0] < s.top()->bounds[0]); //c.a<t.a

			const bool s0_ca_le_tb = s.empty() || (side == 0 && cur_interval->bounds[0] < s.top()->bounds[1]); //c.a < t.b
			const bool s0_tb_le_cb = s.empty() || (side == 0 && s.top()->bounds[1] < cur_interval->bounds[1]); // t.b < c.b

			if (!s.empty() && ((s1_ta_le_cb && s1_ca_le_ta) || ( s0_ca_le_tb && s0_tb_le_cb))) {
				//std::cout << cur_interval << " overlaps " << s.top() << " " << side << std::endl;
				g.new_edge(cur_interval->represented_by[side], s.top()->represented_by[side]); //top becomes parent
			}
			//std::cout << "\tpush " << cur_interval << std::endl;
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
			//std::cout << "glueing equivalence class" <<std::endl;
			slist<node> nodes_to_merge;
			interval<A>* i;
			forall(i, *a_class) {
				assert(i->represented_by[1]==NULL);
				//std::cout << "\t" << i->bounds[0] << " " << i->bounds[1] << " " << i->cont << std::endl;
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
