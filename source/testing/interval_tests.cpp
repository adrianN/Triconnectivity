#include "InterlockingIntervals.hpp"
#include "LEDA/core/slist.h"

void interval_test() {
	interval<int>* intervals[8];
	unsigned int b[8][2] = {
			{0,3},
			{1,3},
			{2,6},
			{2,4},
			{5,8},
			{6,9},
			{7,8},
			{20,25}
	};

	slist<interval<int>* > asc;
	slist<interval<int>* > dsc;
	slist<slist<interval<int>*>* > equiv;


	for(unsigned int i = 0; i<8; i++) {
		interval<int>* it =  new interval<int>(b[i],i);
		intervals[i] = it;
		asc.append(it);
	}

	slist<interval<int>*> one;
	one.append(intervals[0]);
	one.append(intervals[1]);
	one.append(intervals[7]);
	equiv.append(&one);

	dsc.append(intervals[7]);
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
	std::cout <<  res << std::endl;
}
