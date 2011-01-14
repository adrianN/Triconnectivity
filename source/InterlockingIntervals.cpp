#include "LEDA/core/stack.h"
#include "LEDA/core/list.h"
#include <vector>

using namespace leda;
using namespace std;

template <type A> class interval{
public:
	A& cont;
	unsigned int bounds[2];
	interval<A>* parent[2];
	bool leaf_in_F[2];
	list<interval<A>* > children[2];
	interval(unsigned int b[], A& content) : cont(content) {
		bounds[0] = b[0];
		bounds[1] = b[1];
		leaf_in_F[0] = true;
		leaf_in_F[1] = true;
	}
};



// input sorted by [0] asc.
vector<interval<A>* > connect_forest(slist<interval<A> >& input_list, unsigned int which) {

	vector<interval<A>* > input_array(input.size());
	unsigned int i=0;
	for(slist<interval<A> >::item cur_item = input.first(); cur_item!=NULL; cur_item = input.next_item(cur_item)) {
		input_array[i++] = &input_list.contents(cur_item);
	}

	//solve ANLV problem and build forest

	stack<interval<A>* > s;
	for(unsigned int i=input_array.size(); i>=0; i--) {
		interval<A>* cur_interval = input_array[i];

		while(!s.empty() && s.top()->bounds[i] <= cur_interal->bounds[1]) { // < ?
			s.pop();
		}

		if (s.empty()) {
			//no larger value, no parent?
			cur_interval->parent[which] = NULL; //is the root a leaf if no other nodes in the tree
		} else {
			cur_interval->parent[which] = s.top();
			leaf_in_F[which] = false;
			s.top()->children[which].append(cur_interval);
		}

		s.push(cur_interval);
	}
}

void transform_to_rooted_stars(const vector<interval<A>* >& input, unsigned int which) {
	list<interval<A>* > not_done;
	for(unsigned int i = 0; i<input.size(); i++) {
		if (input[i].parent[which]!=NULL)
			not_done.append(input[i]);
	}

	//perform point-to-root pointer jumping
	while(!not_done.empty()) {
		list<interval<A>* > item;
		forall(item, not_done) {
			interval<A>* cur_int = not_done.contents(item);
			assert(cur_int->parent[which] != NULL);

			cur_int->parent[which] = cur_int->parent[which]->parent[which];

			if (cur_int->parent[which]->parent[which] == NULL) {
				not_done.erase(item);
			}
		}
	}
}

// returns vector of roots
vector<interval<A>* > glue_trees(const vector<interval<A>* >& input) {
	for(unsigned int i = 0; i< input.size(); i++) {
		interval<A>* cur_intervall = input[i];

		if (cur_intervall->leaf_in_F[0] && cur_intervall->parent[1] == NULL) {
			interval<A>* the_root_in_F = cur_interval->parent[0];
			interval<A>* child_in_F_prime;
			forall(child_in_F_prime, cur_intervall->children[1]) {
				child_in_F_prime->parent[0] = the_root_in_F;
			}

		}
	}

	vector<interval<A>* > roots;
	for(unsigned int i = 0; i< input.size(); i++) {
		if (input[i]->parent[0] == NULL)
			roots.push_back(input[i]);
	}

	return roots;
}

