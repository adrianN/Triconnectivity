/*
 * testing.hpp
 *
 *  Created on: Dec 5, 2010
 *      Author: adrian
 */

#ifndef TESTING_HPP_
#define TESTING_HPP_

#include "LEDA/graph/ugraph.h"
#include <istream>
#include <memory>

/* test not-3-connected graphs */
void test_separation_pairs(std::istream&);
/* test 3-conn graphs */
void plantri_test(std::istream&);
void plantri_schmidt_test(std::istream&);



std::auto_ptr<leda::ugraph> sample_graph(void);
std::auto_ptr<leda::ugraph> schmidt_sample_graph(void);

std::auto_ptr<leda::ugraph> k4(void);
std::auto_ptr<leda::ugraph> test_eins(void);
std::auto_ptr<leda::ugraph> test_zwei(void);
std::auto_ptr<leda::ugraph> test_drei(void);
std::auto_ptr<leda::ugraph> test_four(void);
std::auto_ptr<leda::ugraph> test_five(void);

#endif /* TESTING_HPP_ */
