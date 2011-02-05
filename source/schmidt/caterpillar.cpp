#include "caterpillar.hpp"

using namespace leda;

caterpillar::caterpillar(void) : parent(NULL) {} //super?


caterpillar::item caterpillar::append(chain* const & element) {
	assert(element!=NULL);
	assert(element->type == two_b);
	if (parent==NULL || element->get_parent()->number<parent->number)
		parent=element->get_parent();
	return slist<chain* const>::append(element);
}
