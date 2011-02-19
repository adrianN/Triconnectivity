#include <memory>
#include "LEDA/graph/ugraph.h"

using namespace std;
using namespace leda;



auto_ptr<ugraph> test_fc_two_children(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,-1,-1},
    		{2,4,-1},
    		{3,0,-1},
    		{1,0,-1},
    		{5,0,-1},
    		{1,0,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
//    {   edge_array<int> edges(*g);
//        edge e;
//        forall_edges(e,*g) {
//            edges[e] = target(e)->id();
//        }
//        g->bucket_sort_edges(edges);
//    }
    return g;
}

auto_ptr<ugraph> test_root_two_children(void) {
	const unsigned int num_nodes = 9;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,4,-1},
    		{2,-1,-1},
    		{3,0,-1},
    		{8,1,-1},
    		{5,-1,-1},
    		{6,-1,-1},
    		{7,4,-1},
    		{5,0,-1},
    		{2,1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
//    {   edge_array<int> edges(*g);
//        edge e;
//        forall_edges(e,*g) {
//            edges[e] = target(e)->id();
//        }
//        g->bucket_sort_edges(edges);
//    }
    return g;
}


auto_ptr<ugraph> test_root_on_one_cycle(void) {
	const unsigned int num_nodes = 9;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,4,-1},
    		{2,-1,-1},
    		{3,0,-1},
    		{8,1,-1},
    		{5,-1,-1},
    		{6,-1,-1},
    		{7,4,-1},
    		{5,4,-1},
    		{2,1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
//    {   edge_array<int> edges(*g);
//        edge e;
//        forall_edges(e,*g) {
//            edges[e] = target(e)->id();
//        }
//        g->bucket_sort_edges(edges);
//    }
    return g;
}

auto_ptr<ugraph> sample_graph(void) {
    auto_ptr<ugraph> g(new ugraph());
    node nodes[13];
    for(unsigned int i=0; i<13; i++) {
        nodes[i] = g->new_node();
    }
    int edges[13][6] = {
        {1,3,7,11,12,-1},
        {12,2,-1,-1,-1,-1},
        {12,3,-1,-1,-1,-1},
        {4,5,6,-1,-1,-1},
        {5,6,7,-1,-1,-1},
        {6,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {10,8,11,-1,-1,-1},
        {9,10,11,-1,-1,-1},
        {10,11,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<13; i++) {
        for(unsigned int j=0; j<6; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

std::auto_ptr<leda::ugraph> schmidt_sample_graph(void) {
	const unsigned int num_nodes = 17;
	const unsigned int max_degree = 7;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{14,16,15,9,3,1,-1}, //v0
    		{6,3,2,-1,-1,-1,-1},
    		{7,5,4,3,-1,-1,-1},
    		{-1,-1,-1,-1,-1,-1,-1},
    		{9,7,6,5,-1,-1,-1}, //v4
    		{6,-1,-1,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1,-1,-1},
    		{8,-1,-1,-1,-1,-1,-1},
    		{16,13,9,-1,-1,-1,-1},
    		{13,10,16,-1,-1,-1,-1}, //v9
    		{15,11,-1,-1,-1,-1,-1},
    		{15,14,12,-1,-1,-1,-1},
    		{13,14,-1,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1,-1,-1}, //v14
    		{-1,-1,-1,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1,-1,-1}

    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }

    return g;

}


auto_ptr<ugraph> k4(void) {
	const unsigned int num_nodes = 4;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,3,2},
    		{2,3,-1},
    		{3,-1,-1},
    		{-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}
/* K4 - one edge */
auto_ptr<ugraph> test_eins(void) {
	const unsigned int num_nodes = 4;
	const unsigned int max_degree = 3;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,3,-1},
    		{2,3,-1},
    		{3,-1,-1},
    		{-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

auto_ptr<ugraph> test_zwei(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 5;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,2,3,-1,-1},
    		{2,3,4,5,-1},
    		{3,4,5,-1,-1},
    		{-1,-1,-1,-1,-1},
    		{5,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

auto_ptr<ugraph> test_drei(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 5;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,2,3,4,-1},
    		{2,3,4,5,-1},
    		{3,4,5,-1,-1},
    		{-1,-1,-1,-1,-1},
    		{5,-1,-1,-1,-1},
    		{-1,-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

// not triconnected
auto_ptr<ugraph> test_four(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 4;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,4,5,-1},
    		{2,3,5,-1},
    		{3,4,-1,-1},
    		{4,-1,-1,-1},
    		{5,-1,-1,-1},
    		{-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}

//triconnected
auto_ptr<ugraph> test_five(void) {
	const unsigned int num_nodes = 6;
	const unsigned int max_degree = 4;
    auto_ptr<ugraph> g(new ugraph());
    node nodes[num_nodes];
    for(unsigned int i=0; i<num_nodes; i++) {
        nodes[i] = g->new_node();
    }
    int edges[num_nodes][max_degree] = {
    		{1,3,4,5},
    		{2,3,4,5},
    		{3,4,-1,-1},
    		{4,5,-1,-1},
    		{5,-1,-1,-1},
    		{-1,-1,-1,-1}
    };

    for(unsigned int i=0; i<num_nodes; i++) {
        for(unsigned int j=0; j<max_degree; j++) {
            if (edges[i][j]>=0) {
                g->new_edge(nodes[i],nodes[edges[i][j]]);
            }
        }
    }
    {   edge_array<int> edges(*g);
        edge e;
        forall_edges(e,*g) {
            edges[e] = target(e)->id();
        }
        g->bucket_sort_edges(edges);
    }
    return g;
}
