/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include <LEDA/core/slist.h>
#include "dfs.hpp"
#include "triconnectivity.hpp"
#include "utilities.hpp"
#include "testing.hpp"
#include "InterlockingIntervals.hpp"
#include <iostream>

using namespace leda;
using namespace std;

void to_file(const ugraph& g, const char* name) {
    std::ostringstream s;
    s << "./graph" << name << ".dot";
    std::fstream f(s.str().c_str(),std::ios_base::out);
    simple_to_dot(g, f);
}






int main(int argc, char* argv[]) {

	interval<int>* intervals[7];
	unsigned int b[7][2] = {
			{0,3},
			{1,3},
			{2,6},
			{2,4},
			{5,8},
			{6,9},
			{7,8}
	};

	slist<interval<int>* > asc;
	slist<interval<int>* > dsc;
	slist<slist<interval<int>*>* > equiv;

	for(unsigned int i = 0; i<7; i++) {
		interval<int>* it =  new interval<int>(b[i],i);
		intervals[i] = it;
		asc.append(it);
	}

	dsc.append(intervals[5]);
	dsc.append(intervals[6]);
	dsc.append(intervals[4]);
	dsc.append(intervals[2]);
	dsc.append(intervals[3]);
	dsc.append(intervals[1]);
	dsc.append(intervals[0]);

	interval<int>* start = intervals[0];
	std::vector<interval<int>*> out;
	bool res = Order<int>::compute_order(asc,dsc,equiv,start,out);
	std::cout << res << std::endl;
    return 0;
}
