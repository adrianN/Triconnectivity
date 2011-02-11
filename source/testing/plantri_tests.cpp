#include "utilities.hpp"
#include "triconnectivity.hpp"
#include "LEDA/graph/ugraph.h"
#include <iostream>
#include <fstream>


using namespace std;
using namespace leda;

void plantri_test(std::istream& cin) {
    ugraph g;
    node s1,s2;
    unsigned int graph_number = 0;
    while(!cin.eof()) {

    	graph_number++;
    	cin >> g;
    	if (!hopcroft_tarjan_is_triconnected_nc(g,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << ".dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(g,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(g,tri);
				tri.close();
			}
			//exit(-1);
    	}
    	if (graph_number % 10000 == 0)
    		cout << '.';
    	if (graph_number % 100000 == 0)
    		cout << '\n';
    	assert(g.number_of_nodes() > 0);
    }

    cout << graph_number << endl;
}

void plantri_schmidt_test(std::istream& cin) {
    ugraph g;
    unsigned int graph_number = 0;
    while(!cin.eof()) {
    	node s1,s2;
    	graph_number++;
    	cin >> g;
//    	std::cout << "************\n\t" << graph_number << "\n********\n" << std::endl;
    	if (!schmidt_is_triconnected(g,s1,s2)) {
    		cout << "drama "<< graph_number << " found separation pair " << s1->id() << " " << s2->id() << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << ".dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(g,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(g,tri);
				tri.close();
			}
			//exit(-1);
    	}
    	if (graph_number % 2500 == 0){
    		cout << '.';
			cout.flush();
    	}
    	if (graph_number % 25000 == 0)
    		cout << " " << graph_number << std::endl;
    	assert(g.number_of_nodes() > 0);
    }

    cout << graph_number << endl;
}


void test_separation_pairs(std::istream& cin) {

    ugraph one,two;
    node s1=NULL,s2=NULL;
    node found1=NULL, found2=NULL;
    unsigned int graph_number = 0;

    while(!cin.eof()) {
    	s1=NULL;
    	s2=NULL;
    	found1=NULL;
    	found2=NULL;
    	graph_number++;
    	cin >> one;
    	if (cin.eof()) return;
    	cin >> two;


    	if (!hopcroft_tarjan_is_triconnected_nc(one,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << "_a.dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(one,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(one,tri);
				tri.close();
			}
			exit(-1);
    	}

    	if (!hopcroft_tarjan_is_triconnected_nc(two,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << "_b.dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(two,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(two,tri);
				tri.close();
			}
			exit(-1);
    	}



    	glue_graphs(one,two,s1,s2);
    	const int id1 = s1->id(), id2 = s2->id();

    	const bool triconnected = hopcroft_tarjan_is_triconnected_nc(one,found1,found2);
    	bool correct_pair = false;
    	if (!triconnected)
    		correct_pair = !((id1==found1->id() && id2==found2->id()) || (id1 == found2->id() && id2 == found1->id()));

    	if ( triconnected || correct_pair) {
    		cout << "drama "<< graph_number << " " << triconnected << " " << correct_pair << endl;
    		if (found1!=NULL && found2!=NULL)
				cout << " found " << found1->id() << " " << found2->id()  << endl;
    		else
    			cout << "\tno pair found" << endl;
    		cout << "\tCorrect pair " << id1 << " " << id2 << endl;
			std::ostringstream s;
			s << "./graph" << graph_number << ".dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(one,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(one,tri);
				tri.close();
			}
			exit(-1);
    	}

    	if (graph_number % 5000 == 0) {
    		cout << '.';
    		cout.flush();
    	}
    	if (graph_number % 50000 == 0)
    		cout << '\n';
    }


}

void schmidt_test_separation_pairs(std::istream& cin) {

    ugraph one,two;
    node s1=NULL,s2=NULL;
    node found1=NULL, found2=NULL;
    unsigned int graph_number = 0;

    while(!cin.eof()) {
    	s1=NULL;
    	s2=NULL;
    	found1=NULL;
    	found2=NULL;
    	graph_number++;
    	cin >> one;
    	if (cin.eof()) return;
    	cin >> two;

    	std::cout << "*****\n\t"<<graph_number<<"\n*******"<<std::endl;

    	if (!hopcroft_tarjan_is_triconnected_nc(one,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << "_a.dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(one,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(one,tri);
				tri.close();
			}
			exit(-1);
    	}

    	if (!hopcroft_tarjan_is_triconnected_nc(two,s1,s2)) {
    		cout << "drama "<< graph_number << endl;
			std::ostringstream s;
			s << "./sep_pair_graph" << graph_number << "_b.dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(two,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(two,tri);
				tri.close();
			}
			exit(-1);
    	}



    	glue_graphs(one,two,s1,s2);
    	const int id1 = s1->id(), id2=s2->id();

    	const bool triconnected = schmidt_is_triconnected(one,found1,found2);
    	bool correct_pair = false;
    	if (!triconnected)
    		correct_pair = !((id1==found1->id() && id2==found2->id()) || (id1 == found2->id() && id2 == found1->id()));

    	if ( triconnected || correct_pair) {
    		cout << "drama "<< graph_number << " " << triconnected << " " << correct_pair << endl;
    		if (found1!=NULL && found2!=NULL)
				cout << " found " << found1->id() << " " << found2->id()  << endl;
    		else
    			cout << "\tno pair found" << endl;
    		cout << "\tCorrect pair " << id1 << " " << id2 << endl;
    		hopcroft_tarjan_is_triconnected_nc(one,s1,s2);
    		cout << "\tHT says " << s1->id() << " "<< s2->id() << endl;
			std::ostringstream s;
			s << "./graph" << graph_number << ".dot";
			{	std::fstream dot(s.str().c_str(), std::ios::out);
				simple_to_dot(one,dot);
				dot.close();
			}
			s.seekp((long)s.tellp()-4);
			s << ".tri";
			{	std::fstream tri(s.str().c_str(),std::ios::out);
				write_planar_code(one,tri);
				tri.close();
			}
			exit(-1);
    	}

    	if (graph_number % 5000 == 0) {
    		cout << '.';
    		cout.flush();
    	}
    	if (graph_number % 50000 == 0)
    		cout << '\n';
    }


}
