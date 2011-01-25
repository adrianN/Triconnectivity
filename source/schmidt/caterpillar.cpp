#include "caterpillar.hpp"

using namespace leda;

caterpillar::caterpillar(void) : parent(NULL) {} //super?
caterpillar::item caterpillar::append(chain* const & element) {
	assert(element!=NULL);
	if (parent==NULL || element->number<parent->number)
		parent=element;
	return slist<chain* const>::append(element);
}
