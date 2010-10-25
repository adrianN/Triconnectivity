/*
 * dictionary.h
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

struct item {
	void const * key ;
	void const * value ;
};

typedef item* dic_item;

template <typename K, typename I> class dictionary {

public:

	dictionary() {};
	virtual ~dictionary() {
		for(unsigned int i=0; i<pairs.size(); i++) {
			delete pairs[i];
		}
	};

	dic_item insert(const K& key, const I& value) {
		//std::cout << "insert " << key << " " << value << std::endl;
		dic_item contained = lookup(key);
		if (contained!=NULL) {
			contained->value=&value;
			return contained;
		} else {
			item* it = new item;
			it->key=&key;
			it->value=&value;
			pairs.push_back(it);
			std::cout << "insert " << it->key << " " << it->value << std::endl;
			return it;
		}
	};
	const K& key(const dic_item& it) const {
		return *static_cast<K const *>(it->key);
	};

	const I& inf(const dic_item& it) const {
		return *static_cast<I const *>(it->value);
	};

	item* lookup(const K& key) {
		int i = pos(key);
		if (i>=0) return pairs[i];
		return NULL;
	}

	void del(const K& key) {
		pairs.erase(pairs.begin()+pos(key));
	}

	int size() const {
		return pairs.size();
	}

	const dic_item& operator[](int i) const {
		return pairs[i];
	}

	std::string toString() const {
		std::ostringstream s;
		for(int i =0; i<size(); i++) {
			dic_item it = pairs[i];
			s << "(" << it->key << "," << it->value << "), ";
		}
		return s.str();
	}
private:
	std::vector<item*> pairs;

	int pos(const K& key) const {
		for(unsigned int i=0; i<pairs.size(); i++) {
			if ((pairs[i]->key)==&key) {
				return i;
			}
		}
		return -1;
	};

};

#endif /* DICTIONARY_H_ */
