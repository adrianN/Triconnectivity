#include "triconnectivity.hpp"
#include <ostream>

using namespace leda;
using std::endl;

void to_dot(const ugraph& g, std::ostream& out) {

	out << "digraph G {" << endl;
	node n;
	forall_nodes(n,g) {

		out << "\tnode" << n->id() << " [label = \"" << n->id() << "\"]" << endl;

		node u;
		unsigned int i=0;
		forall_adj_nodes(u,n) {
			out << "\tnode" << n->id() << " -> " << "node" << u->id() << " [label = \"" << ++i << "\"]" << endl;
		}
	}

	out << "}" << endl;

}
