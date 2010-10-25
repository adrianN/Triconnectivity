/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include "dfs.hpp"
#include "triconnectivity.hpp"
#include "utilities.hpp"
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

using namespace leda;
using std::cout;
using std::endl;
using std::ostream;

int main(void) {
	ugraph g;
	const double p = 0.5;
	const unsigned int n = 10;
	const unsigned int m = (int)(((double)(n*(n-1))*p)/2.0);
	cout << "n, m " << n << ", "<<m<<endl;
	float total = used_time();
	for(int i=0; i<1; i++) {
		random_graph(g,n,m,false,true,false);

		std::ostringstream s;
		s << "./graph" << i << ".dot";
		std::fstream f(s.str().c_str(),std::ios_base::out);
		to_dot(g, f);
		cout << "==" <<i << "==" <<endl;
		cout << is_connected(g) << endl;
		cout << naive_is_triconnected(g) << endl;


	}
	cout << "total " << used_time()-total << endl;
	return 0;
}
