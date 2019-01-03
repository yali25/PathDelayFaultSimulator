/*
   This file contains several typdes which are used in the program. 
   NOTE: the def. typedef boost::dynamic_bitset<unsigned long long> boost_bitset;
   defines a bitset which uses a 64bit unsigned integer to store the bits.
   This bitset is used through out the calculation of the sequential version of the algorithm.
*/
#ifndef HASH_MAPS
#define HASH_MAPS

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <string>

typedef boost::dynamic_bitset<unsigned long long> boost_bitset;

// this hashmap is used to represent a node in the graph, the first string is the name of the node from the netlist
// the pair saves the gate where it appeared first and a unique inetger
typedef std::unordered_map<std::string, std::pair<std::string, std::size_t>>my_map;
typedef std::unordered_map<std::string,std::pair<std::pair<std::string,std::string>, std::pair<std::string,std::string>>> dff_map;
typedef std::unordered_map<std::string, boost_bitset> bitset_map;




#endif 