// This file contains the typedef for the graph
// for more details check the Boost Graph Library doc
#ifndef GRAPH_STRUCTURE
#define GRAPH_STRUCTURE

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;

#endif
