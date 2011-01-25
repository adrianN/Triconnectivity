#include "certificate.hpp"

certificate::certificate(ugraph& graph, unsigned int k, cert_data* c) : the_graph(graph), connectivity(k), content(c) {assert(k<=3);}
bool certificate::verify() {
	return true;
}
