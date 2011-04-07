/*
 * dictionary.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: adrian
 */

#include <LEDA/graph/ugraph.h>
#include <LEDA/graph/node_array.h> //never include before graph.h
#include <LEDA/core/slist.h>
#include <LEDA/core/list.h>
#include <LEDA/geo/d3_rat_point.h>
#include <LEDA/geo/d3_hull.h>
#include <LEDA/core/random_source.h>
#include "dfs.hpp"
#include "triconnectivity.hpp"
#include "utilities.hpp"
#include "testing.hpp"
#include "InterlockingIntervals.hpp"
#include "LEDA/system/timer.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <LEDA/geo/point_set.h>
#include <LEDA/geo/point.h>
#include <LEDA/geo/point_dictionary.h>

using namespace leda;
using namespace std;

void check_intervals() {
	unsigned int (*fc)(interval<int>* const &) = &interval<int>::first_component;
	unsigned int (*sc)(interval<int>* const &) = &interval<int>::second_component;

	interval<int> start(0,2,0);

	slist<interval<int>*> intervals;
	intervals.append(&start);
	intervals.append(new interval<int>(1,3,1));
	intervals.append(new interval<int>(1,3,2));
	intervals.append(new interval<int>(1,3,3));

	slist<interval<int>*> intervals_asc = bucket_sort<interval<int>*,asc, dsc>(intervals,fc,sc,0,3);
	slist<interval<int>*> intervals_dsc = bucket_sort<interval<int>*,dsc, asc>(intervals,sc,fc,0,3);
	interval<int>* n;
	forall(n, intervals_asc) {
		std::cout << n << " ";
	}
	std::cout << std::endl;
	forall(n,intervals_dsc) {
		std::cout << n << " ";
	}
	std::cout << std::endl;
	slist<slist<interval<int>*>* > eq;
	std::vector<interval<int>*> out;
	Order<int>::compute_order(intervals_asc, intervals_dsc ,eq , &start, out);
	for(unsigned int i =0; i<out.size(); i++)
		std::cout << out[i] << std::endl;

}

void to_file(const ugraph& g, const char* name) {
    std::ostringstream s;
    s << "./graph" << name << ".dot";
    std::fstream f(s.str().c_str(),std::ios_base::out);
    simple_to_dot(g, f);
}

bool check_graph(ugraph u, timer& htt, timer& st, bool both) {

	//simple_to_dot(u,f);
	//std::cout << "n " << u.number_of_nodes() << " m " << u.number_of_edges() << std::endl;
	node s1=NULL, s2=NULL;
	node s3=NULL, s4=NULL;
	bool s = false, ht = false;
	if (both) {
		st.start();
		s =  schmidt_is_triconnected(u,s1,s2);
		st.stop();
	}
	htt.start();
	ht = hopcroft_tarjan_is_triconnected(u,s3,s4);
	htt.stop();
	if (both) {
		bool all = s == ht;
		if (!all) {
			bool naive = naive_is_triconnected(u);

			to_file(u,"the_graph");
			fstream f("./the_graph.tri", ios::out);
			write_planar_code(u,f);

			std::cout << "ht " << ht << " s " << s << " naive " << naive << std::endl;
			std::cout << "Sep ";
			if (s1!=NULL) {
				u.hide_node(s1);
				std::cout << s1->id() << " ";
			} else {
				std::cout << "null ";
			}
			if (s2!=NULL) {
				u.hide_node(s2);
				std::cout << s2->id() << " ";
			} else {
				std::cout << "null ";
			}
			std::cout << std::endl;
			to_file(u,"hidden");
			assert(false);
		}
	}

//	std::cout << "ht " << ht << " s " << s << std::endl;
	return ht || s;
}

bool check_random_graph(int n, int m, timer& ht, timer& s, bool both) {
	graph g;
	random_simple_undirected_graph(g,n,m);
	ugraph u(g);
	return check_graph(u,ht,s, both);
}

void f() {
	//	fstream f;
	//
	//	if (argc>1)
	//		{
	//			cout << "reading file " << argv[1] << std::endl;
	//
	//			f.open(argv[1], std::ios::in);
	//		}
	//	else {
	//		cout << "reading file ./9.txt" << std::endl;
	//		f.open("./9.txt", std::ios::in);
	//	}
	//
	////	ugraph g;
	//////	g.read("./the_graph.leda");
	////	f >> g;
	////	g.permute_edges();
	////	node s3 = NULL,s4 = NULL;
	////	hopcroft_tarjan_is_triconnected_nc(g,s3,s4);
	////	node n;
	////	forall_nodes(n,g) {
	////		node s1 = NULL, s2 = NULL;
	////
	////		std::cout << "Starting node " << n->id() << std::endl;
	////	schmidt_is_triconnected(g,s1,s2,n);
	////		if (s1!=NULL && s2 != NULL) {
	////			std::cout <<s1->id() << " " << s2->id() << std::endl;
	////			std::cout << s3->id() << " " << s4->id() << std::endl;
	////			if (!((s1->id() == s3->id() && s2->id() == s4->id()) || (s1->id() == s4->id() && s2->id() == s3->id()))) {
	////				exit(-1);
	////			}
	////		} else {
	////			std::cout << "none found" << std::endl;
	////
	////			std::cout << (s3!=NULL ? s3->id() : -1) << " " << (s4!=NULL? s4->id():-1) << std::endl;
	////			exit(-1);
	////		}
	////	}
	//
	////	plantri_test(f);
	//	plantri_schmidt_test(f);
}

void check_planar_code(const char* filename) {
	fstream f(filename, ios::in);
	ugraph u;
	f >> u;
	node s1,s2;
	bool s =  schmidt_is_triconnected(u,s1,s2);
	bool ht = hopcroft_tarjan_is_triconnected(u,s1,s2);
	if (ht!=s) {
		std::cout << "ht " << ht << " s " << s << " naive " << naive_is_triconnected(u) << std::endl;
	} else {
		std::cout << "a-okay" << std::endl;
	}
}

float prob(int m, int num, int n, bool both) {
//	const int edges = (n*(n-1))/2;
	int count=0;
//	const float p = f/((float)n);
//	const int m = (int)((float)edges*p);
	for (int i = 0; i<num; i++) {
		timer ht, s;
		bool c = check_random_graph(n,m,ht,s,both);
		if (c) {
			count++;
		}
		std::cout << n << " " << num << " " << m << " "  << (m+n) << " " << c << " " << ht.elapsed_time() << " " << s.elapsed_time()  <<std::endl;

	}

	return (float)count/num;
}

//void b_search(float l, float r, int n) {
//	float m = l+(r-l)/2;
//
//	if ((r-l)<0.2) {
//		prob(m,n/100,n,true);
//		return;
//	}
//	float p = prob(m,std::max(50,(int)(30/(r-l))),n,false);
//	if (p>0.5)
//		b_search(l,m,n);
//	else
//		b_search(m,r,n);
//}

void graph_statistics(int n) {
	float max = 0.5*log(n)+1.5*log(log(n))-0.7;
	prob(max*n,10,n,true);
//	float f=0;
//	while(f<0.5) {
//		max*=2;
//		f = prob(max,50,n,false);
//	}
//	b_search(max/1.5,max,n);
//	const int edges = (n*(n-1))/2;

//	for (float f = 2*log(n); f < n/2; f+=n/200) {
//		timer ht,s;
//		int count=0;
//		const float p = f/((float)n);
//		const int m = (int)((float)edges*p);
//		std::cout  << n << " " << f << " ";
//		for (int i = 0; i<num; i++) {
//
//			if (check_random_graph(n,m,ht,s)) {
//				count++;
//			}
//		}
//		std::cout << (float)count/num << " " << (ht.elapsed_time()/(float)num) << " " << (s.elapsed_time()/(float)num) << " " << (s.elapsed_time()/ht.elapsed_time()) <<std::endl;
//	}
}

void add_edges() {
	ugraph g;
	std::cout << "n m d\n";
	for(int num_nodes = 1000; num_nodes < 1000000; num_nodes+=1000) {
		int m = num_nodes;
		for (int k=0; k<200; k++) {
			random_simple_undirected_graph(g,num_nodes,m);
			std::vector<node> nodes;
			node n;
			forall_nodes(n,g) {
				nodes.push_back(n);
			}
			random_source r(0,nodes.size()-1);
			bool added=false;
			do {
				added=false;
				node n;
				forall_nodes(n,g) {
					if (g.degree(n) <3) {
						added=true;
						node u = n;
						while(u==n || are_connected(n,u)) {
							int un;
							r >> un;
							u = nodes[un];
						}
						g.new_edge(n,u);
						break;
					}
				}

			} while(added);
			forall_nodes(n,g) {
				if (g.degree(n) < 3) {
					std::cout << "DRAMA";
					exit(0);
				}
			}
			int d = 0;
			node s1,s2;
			while(!hopcroft_tarjan_is_triconnected_nc(g,s1,s2)) {
				int un=0,vn=0;
				while(un==vn) {
				r >> un;
				r >> vn;
				}
				g.new_edge(nodes[un],nodes[vn]);
				d++;
			}
			std::cout << g.number_of_nodes() << " " << g.number_of_edges() << " " << d << std::endl;
		}
	}
}


void remove_edges() {
	for(unsigned int num_nodes=4; num_nodes<10; num_nodes+=1) {
		for(int k =0; k < 10 ; k++) {
			ugraph g;
			node s1,s2;
			bool tr=false;
//			int m = num_nodes*(0.5*log(num_nodes)+1.5*log(log(num_nodes)) -0.8);
			int m = ((num_nodes-1)*num_nodes)/2;
			do {
				g.clear();
				random_simple_undirected_graph(g,num_nodes,m);
				tr = hopcroft_tarjan_is_triconnected_nc(g,s1,s2);
//				std::cout << tr << std::endl;
			} while(!tr);
			to_file(g,"ausgang");
			std::vector<edge> edges;
			edge e;
			forall_edges(e,g) {
				if (g.degree(source(e)) >3 && g.degree(target(e)) > 3) {
					edges.push_back(e);
				}
			}
			std::random_shuffle(edges.begin(), edges.end());

			bool triconnected = true;
			unsigned int i=0;
			int c = 1;
			while(triconnected) {
				edge e=NULL;
				while(i<edges.size() && (g.degree(source(edges[i]))<=3 || g.degree(target(edges[i]))<=3)) {
					i++;
				}
				if (i==edges.size()) {
					edge e;
					forall_edges(e,g) {
						if (!(g.degree(source(e))==3 || g.degree(target(e)) == 3)) {
							std::cout << "badumm";
							exit(0);
						}
					}
					c=0;
					break;
				}
				e = edges[i];
				i++;
				assert(g.degree(source(e)) > 3 && g.degree(target(e)) > 3);
				g.del_edge(e);
				std::cout << "del " << source(e)->id() << " " << target(e) ->id() << std::endl;
//				std::cout << g.number_of_edges() << std::endl;
				triconnected = hopcroft_tarjan_is_triconnected_nc(g,s1,s2);
			}
			if (!triconnected) {
				to_file(g,"uebrig");
				std::cout << "oho" << std::endl;
				exit(0);
			}
			std::cout << g.number_of_nodes() << " "  << g.number_of_edges() << " " << c << " " << triconnected << std::endl;
		}
	}
}

void check_convex_hull(int n) {
	random_source r;
	for(unsigned int i=0; i<50; i++) {
	list<d3_rat_point> points;

	random_d3_rat_points_on_sphere(n,100000, points);

		GRAPH<d3_rat_point,int> h;
		CONVEX_HULL(points,h);



		ugraph g(h);
		make_simple(g);
		timer ht,s;
		std::cout << g.number_of_nodes() << " " << g.number_of_edges() << " ";
		check_graph(g,ht,s,true);
		std::cout << ht.elapsed_time() << " " << s.elapsed_time() << std::endl;
	}

}

void check_dense_graphs(int n) {
	int m = sqrt(n)*n;
	node x,y;
	int num = 10;
	bool both=true;
	int count=0;

	for (int i = 0; i<num; i++) {
		graph g,h;
		random_simple_undirected_graph(g,n,m);
		random_simple_undirected_graph(h,n,m);
		ugraph u(g);
		ugraph v(h);
		glue_graphs(u, v, x, y);

		timer ht, s;
		bool c = check_graph(u,ht,s,both);
		if (c) {
			count++;
		}
		std::cout << u.number_of_nodes() << " " << num << " " << u.number_of_edges() << " "  << (u.number_of_nodes()+u.number_of_edges()) << " " << c << " " << ht.elapsed_time() << " " << s.elapsed_time()  <<std::endl;
	}
}

void density(int m) {
	node x,y;
	int num = 10;
	bool both=true;
	int count=0;
	int n = 20000;

	for (int i = 0; i<num; i++) {
		graph g,h;
		random_simple_undirected_graph(g,n,m);
		random_simple_undirected_graph(h,n,m);
		ugraph u(g);
		ugraph v(h);
		glue_graphs(u, v, x, y);

		timer ht, s;
		bool c = check_graph(u,ht,s,both);
		if (c) {
			count++;
		}
		std::cout << u.number_of_nodes() << " " << num << " " << u.number_of_edges() << " "  << (u.number_of_nodes()+u.number_of_edges()) << " " << c << " " << ht.elapsed_time() << " " << s.elapsed_time()  <<std::endl;
	}
}

void density_one_graph(int m) {
	int num = 7;
	bool both=true;
	int count=0;
	int n = 40000;

	for (int i = 0; i<num; i++) {
		graph g,h;
		random_simple_undirected_graph(g,n,2*m);
		ugraph u(g);

		timer ht, s;
		bool c = check_graph(u,ht,s,both);
		if (c) {
			count++;
		}
		std::cout << u.number_of_nodes() << " " << num << " " << u.number_of_edges() << " "  << (u.number_of_nodes()+u.number_of_edges()) << " " << c << " " << ht.elapsed_time() << " " << s.elapsed_time()  <<std::endl;
	}
}

void geometric_graphs(int n) {
	for(int k=1; k<log(n)+1; k++) {
	for(int i=0; i<200; i++) {
		list<point> points;
		random_points_in_square(n,n,points);
		d2_dictionary<double,double,node> nodes;
		ugraph g;
		{point p;
		forall(p,points) {
			nodes.insert(p.xcoord(),p.ycoord(), g.new_node());
		}
		}
		point_set set(points);
		point p;
		forall(p,points) {
			list<node> crappy_neighbours = set.nearest_neighbors(p,k+1);
			node n;
			crappy_neighbours.pop();
			list<point> neighbours;
			forall(n,crappy_neighbours) {
				neighbours.append(set.pos(n));
			}
			point q;
			forall(q,neighbours) {
				g.new_edge(nodes[nodes.lookup(p.xcoord(),p.ycoord())],nodes[nodes.lookup(q.xcoord(),q.ycoord())]);
			}
		}
		make_simple(g);
		//to_file(g,"blubb");
		node s1=NULL,s2=NULL;
		bool tconn = hopcroft_tarjan_is_triconnected_nc(g,s1,s2);
		int conn = tconn? 3 : (s1!=s2 ? 2 : (s1 !=NULL ? 1 : 0));
		std::cout << g.number_of_nodes() << " " << g.number_of_edges() << " " << k << " " << tconn <<  " " << conn << "\n";
	}
	}
}

int main(int argc, char* argv[]) {
	int n = 100;
	if (argc>1) {
		istringstream i(argv[1]);
		i >> n;
	}

//	geometric_graphs(n);
//	remove_edges();
	add_edges();
//	check_intervals();
//	check_dense_graphs(n);
//	density(n);
//	density_one_graph(n);
//	check_convex_hull(n);
//	graph_statistics(n);
//	timer ht,s;
//	check_random_graph(1000,20000,ht,s);
//	std::cout << ht.elapsed_time() << " " << s.elapsed_time() << std::endl;
return 0;
}
