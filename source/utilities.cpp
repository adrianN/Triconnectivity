#include "triconnectivity.hpp"
#include <ostream>
#include <sstream>

using namespace leda;
using std::endl;

void to_dot(const ugraph& g, std::ostream& output) {
	edge e;
	std::ostringstream out;
	out << "graph G {" << endl;

	forall_edges(e,g) {
		node src = source(e);
		node trg = target(e);

		out << "\tnode";
		g.print_node(src,out);
		out << " -- node";
		g.print_node(trg,out);
		out << endl;
	}

	out << "}" << endl;
	std::string s = out.str();
	// Since dot doesn't like [], we have to replace that. There should be a fancy way to do global replace in strings
	for(unsigned int i=0; i<s.length(); i++) {
		if (s[i]=='[' || s[i] == ']') {
			s[i] = '_';
		}
	}

	output << s;
}
