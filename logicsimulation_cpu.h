/*
	This file contains a simple logic simulation for the CPU with 2 patterns.
	This logic simulation is NOT used by the path classification algorithm.
	This logic simulation is used when doing a stilefile check.
*/

#ifndef LOGICSIM_FILE
#define LOGICSIM_FILE

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include "graph_structure.h"
#include "helper_functions.h"
#include "hash_maps.h"
#include <vector>


void logic_simlation_cpu(const Graph &g,
					     const std::vector<unsigned int> &top_sorted_graph,
						 const std::vector<char> &gate_types_vec,
						 const boost_bitset::size_type    &bppk,
						 std::vector<boost_bitset> &circuit_states1,
						 std::vector<boost_bitset> &circuit_states2
						 )
{
	boost_bitset result1(bppk);
	boost_bitset result2(bppk);

	for (unsigned int topsorted_iterator=0; topsorted_iterator < top_sorted_graph.size(); topsorted_iterator++)
	{
		if (boost::in_degree(topsorted_iterator, g))
		{
			auto inEdges = boost::in_edges(topsorted_iterator, g);  // Get an iterator to the incident edges of the current node

			if (gate_types_vec[topsorted_iterator] == 6) // Gate is an inverter
			{
				result1 = ~circuit_states1[(*inEdges.first).m_source];
				result2 = ~circuit_states2[(*inEdges.first).m_source];
			}
			else // All the other gates with 2 or more inputs
			{
				result1 = circuit_states1[(*inEdges.first).m_source];
				result2 = circuit_states2[(*inEdges.first).m_source];
				++inEdges.first;
				if (gate_types_vec[topsorted_iterator] == 1 || gate_types_vec[topsorted_iterator] == 4) // AND and NAND
				{

					for (; inEdges.first != inEdges.second; ++inEdges.first)
					{
						result1 &= circuit_states1[(*inEdges.first).m_source];
						result2 &= circuit_states2[(*inEdges.first).m_source];
					}

					if (gate_types_vec[topsorted_iterator] == 1) // if NAND then invert the results
					{
						result1.flip();
						result2.flip();
					}

				}
				else if (gate_types_vec[topsorted_iterator] == 5 || gate_types_vec[topsorted_iterator] == 2) // OR and NOR casse
				{

					for (; inEdges.first != inEdges.second; ++inEdges.first)
					{
						result1 |= circuit_states1[(*inEdges.first).m_source];
						result2 |= circuit_states2[(*inEdges.first).m_source];
					}
					if (gate_types_vec[topsorted_iterator] == 2) // if NOR then invert the results
					{
						result1.flip();
						result2.flip();
					}


				}
				else if (gate_types_vec[topsorted_iterator] == 3)  // XOR
				{
					for (; inEdges.first != inEdges.second; ++inEdges.first)
					{
						result1 ^= circuit_states1[(*inEdges.first).m_source];
						result2 ^= circuit_states2[(*inEdges.first).m_source];
					}
				}
			}

			circuit_states1[topsorted_iterator] = std::move(result1);
			circuit_states2[topsorted_iterator] = std::move(result2);
			
		}

	}

}
#endif