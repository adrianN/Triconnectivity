#include "certificate.hpp"
#include "chain_edge_iterator.hpp"


certificate::certificate(ugraph const & graph, const  schmidt_triconnectivity* d) : the_graph(graph), decomposition(d) {}


bool certificate::add_bg_path(list<edge> const & edges) {
	std::cout << "\t****** BG PATH: " ;
	edge e;
	forall(e,edges) {
		std::cout << decomposition->dfi(source(e)) << " -> " << decomposition->dfi(target(e)) << " ";
	}
	std::cout << std::endl;
	return true;
}

bool certificate::add_bg_path(const chain * c) {
	list<edge> l;
	chain_edge_iterator it(c, decomposition);
	for(edge e = it.next(); e!=NULL; e=it.next()) {
		l.append(e);
	}
	add_bg_path(l);
	return true;
}

bool certificate::verify() const {
	return true;
}
