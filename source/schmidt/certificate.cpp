#include "certificate.hpp"

certificate::certificate(ugraph& graph) : the_graph(graph) {}
bool certificate::add_bg_path(std::vector<edge> const & edges) {
	return true;
}

bool certificate::add_bg_path(list<edge> const & edges) {
	return true;
}

bool certificate::add_bg_path(const chain * c) {
	return true;
}

bool certificate::verify() const {
	return true;
}
