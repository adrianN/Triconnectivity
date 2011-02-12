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
#include "not_triconnected_exception.hpp"
#include <vector>
#include <string>
#include <fstream>
using namespace std;
using namespace leda;


//#define PRINT
//#define OUTPUT_DOT

template <typename A> class interval{
public:
	A cont;
	unsigned int bounds[2];
	leda::node represented_by;
	interval(unsigned int l, unsigned int r, A content) : cont(content) {
#ifdef PRINT
		std::cout << "Fresh interval " << l << " " << r  << " " << cont << std::endl;
#endif
		if (l>r)
			std::swap(l,r);
		bounds[0] = l;
		bounds[1] = r;
		represented_by = NULL;
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

	static void compute_order(const slist<interval<A>*>& sorted_asc,
							  const slist<interval<A>*>& sorted_dsc,
							  const slist<slist<interval<A>*>* >& equivalent_intervals,
							  const interval<A>* start,
							  std::vector<interval<A>* >& output) throw(std::pair<unsigned int,unsigned int>)
	{
#ifdef OUTPUT_DOT
		static unsigned char b='0';
#endif
		assert(sorted_asc.size() == sorted_dsc.size());
		assert(sorted_asc.size() > 0);
#ifdef PRINT
		std::cout << "Computing overlap graph " << sorted_asc.size() << std::endl;
#endif
		ugraph g;
		connect_forest(g, sorted_asc, 1);
#ifdef OUTPUT_DOT
		{
			std::string s("a_forest");
			s += b;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}
#endif

		connect_forest(g, sorted_dsc, 0);
#ifdef OUTPUT_DOT
		{
			std::string s("b_forest");
			s += b;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}
#endif

		glue_equivalence_classes(g, equivalent_intervals);
#ifdef OUTPUT_DOT
		{
			std::string s("final_graph");
			s += b++;
			s += ".dot";

			std::fstream f (s.c_str(), std::ios::out);
			to_dot<interval<A>*>(g, labels(g,sorted_asc),f);
			f.close();
		}
#endif

		//mapping from nodes to intervals
		node_array<interval<A>*> belongs_to(g, NULL);
		{
			interval<A>* i;
			forall(i, sorted_asc) {
				belongs_to[i->represented_by] = i;
			}
		}

		node* ordered = new node[g.number_of_nodes()];
		for(unsigned int i=0; i<(unsigned int)g.number_of_nodes(); i++) {
			ordered[i]=NULL;
		}
		const node startnode = start->represented_by;
		bool is_connected = dfs_order(g, startnode, ordered);
		if (!is_connected) {
#ifdef PRINT
			std::cout << "Interlock graph not connected. Main component contains" << std::endl;
#endif
			// we started at the merged nodes of the real vertices. Hence any connected component that is not visited contains the information to find a separation pair.
			//find another connected component
			node_array<bool> visited(g,false);
			for(unsigned int i = 0; i<(unsigned int)g.number_of_nodes() && ordered[i] != NULL; i++) {
					visited[ordered[i]] = true;
#ifdef PRINT
					std::cout << belongs_to[ordered[i]] << " ";
#endif
					ordered[i]=NULL;
			}
			//std::cout << std::endl;
			{	node n;
				forall_nodes(n,g) {
					if (!visited[n]) { //an unvisited connected component
						dfs_order(g,n,ordered);
						break;
					}
				}
			}
			//find the extreme attachments
			unsigned int min = UINT_MAX;
			unsigned int max = 0;
			{
#ifdef PRINT
				std::cout << "One remaining component contains " << std::endl;
#endif
				for(unsigned int i = 0; i<(unsigned int)g.number_of_nodes() && ordered[i] != NULL; i++) {
					interval<A>* it = belongs_to[ordered[i]];
#ifdef PRINT
					std::cout << it << " ";
#endif
					min = std::min(it->bounds[0],min);
					max = std::max(it->bounds[1],max);
				}
#ifdef PRINT
				std::cout <<"\n Min/max " << min << " " << max << std::endl;
#endif
			}
			delete[] ordered;
			throw pair<unsigned int,unsigned int>(min,max);
		}


		for(unsigned int i=0; i<(unsigned int)g.number_of_nodes(); i++) {
			output.push_back(belongs_to[ordered[i]]);
		}
		delete[] ordered;
	}

private:
	static node_array<interval<A>*> labels(const ugraph& g, const slist<interval<A>*> list) {
		node_array<interval<A>*> l(g);
		interval<A>* val;
		forall(val, list) {
			if(val->represented_by != NULL)
				l[val->represented_by] = val;
		}
		return l;
	}

	static void connect_forest(ugraph& g, const slist<interval<A>* >& input_list, unsigned int side) {
		std::vector<interval<A>*> input_array;
		{	// create a node for each interval and copy the intervals into a vector
			interval<A>* val=NULL;
#ifdef PRINT
			std::cout << "side " << side << " intervals: ";
#endif
			forall(val, input_list) {
				assert(val!=NULL);
#ifdef PRINT
				std::cout << val << " ";
#endif
				input_array.push_back(val);
				if (val->represented_by == NULL) {
					assert(side==1);
					val->represented_by = g.new_node();
				}
			}
#ifdef PRINT
			std::cout << std::endl;
#endif
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
				g.new_edge(cur_interval->represented_by, s.top()->represented_by); //top becomes parent
			}
			//std::cout << "\tpush " << cur_interval << std::endl;
			s.push(cur_interval);
		}
	}

	static void glue_equivalence_classes(ugraph& g, const slist<slist<interval<A>*>* > equivalent_intervals) {
		slist<interval<A>* >* a_class;
		forall(a_class,equivalent_intervals) {
			//std::cout << "glueing equivalence class" <<std::endl;
			slist<node> nodes_to_merge;
			interval<A>* i;
			unsigned int min_bound = UINT_MAX;
			unsigned int max_bound = 0;
			forall(i, *a_class) {
				assert(i->represented_by!=NULL);
				//to properly compute the maximal dependent path we need to update the bounds
				min_bound = std::min(min_bound,i->bounds[0]);
				max_bound = std::max(max_bound,i->bounds[1]);
				//std::cout << "\t" << i->bounds[0] << " " << i->bounds[1] << " " << i->cont << std::endl;
				nodes_to_merge.append(i->represented_by);
			}
			node merged_node = merge_nodes(g, nodes_to_merge);
			forall(i,*a_class) {
				i->bounds[0] = min_bound;
				i->bounds[1] = max_bound;
				i->represented_by = merged_node;
			}
		}

	}
};

#endif /* INTERLOCKINGINTERVALS_HPP_ */
