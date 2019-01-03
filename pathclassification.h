/*
	This file contains the path classifiaction algorithm. It follows exactly the description given in
	the Master's thesis and in the two papers (FTC, ITC).
	This file contains three implementenations of the algorithm but only the function "path_classification_with_logic" should
	be used. The logic simualtion is already included in this function. The version "path_classification" does not include a logic
	simulation. 
	The other two versions were developed in order to debug the OpenCL version. They should not be used.
*/


#ifndef PATH_CLASS
#define PATH_CLASS


#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include "graph_structure.h"
#include <vector>


boost_bitset isRobustAndNand(const boost_bitset &i1v1,
	const boost_bitset &i1v2,
	const boost_bitset &i2v1,
	const boost_bitset &i2v2)
{

	//return ((~i1v1) & i1v2 & i2v2) | (i1v1 & (~i1v2) & i2v1 & i2v2);
	return (((~i1v1) & i1v2) | (i1v1 & (~i1v2) & i2v1)) & i2v2;
}

boost_bitset isRobustOrNor(const boost_bitset &i1v1,
								  const boost_bitset &i1v2,
								  const boost_bitset &i2v1,
								  const boost_bitset &i2v2)
{

	//return (i1v1 & (~i1v2) & (~i2v2)) | ((~i1v1) & i1v2 & (~i2v1) & (~i2v2));
	return ((i1v1 & (~i1v2)) | ((~i1v1) & i1v2 & (~i2v1))) & (~i2v2);
}

boost_bitset isRobustOrNonRobustAndNand(const boost_bitset &i1v1,
											 const boost_bitset &i1v2,
											 const boost_bitset &i2v2)
{
	return (i1v1 ^i1v2)& i2v2;
}

boost_bitset isRobustOrNonRobustOrNor(const boost_bitset &i1v1,
										         const boost_bitset &i1v2,
										         const boost_bitset &i2v2)
{
	return (i1v1 ^i1v2)&(~i2v2);
}


// the data structure is_robust_or_non_robust_result saves the Robust or Non Robust information, the unsigned int value is the on input
void path_classification(const Graph                           &g, 
						 const Graph                           &transponse_graph,
						 const std::vector<boost_bitset> &circuit_states1,
						 const std::vector<boost_bitset> &circuit_states2,
						 const std::vector<char>                   &gate_types_vec,
						 const unsigned int                        &counter,
						 const boost_bitset::size_type    &bppk,
						 const std::vector<std::string>             &vertex_names,
						 std::vector<std::vector<std::pair<unsigned int,boost_bitset>>> &is_robust_or_non_robust_result,
						 std::vector<std::vector<std::pair<unsigned int,boost_bitset>>> &is_robust_result,
						 std::vector<boost_bitset> &is_robust_path_to_po,
						 std::vector<std::vector<std::pair<unsigned int,boost_bitset>>> &belongs_to_robust_path_final_result)
{
	std::cout << "Path Classification starts..." << std::endl;
	is_robust_or_non_robust_result.reserve(counter);
	is_robust_result.reserve(counter);
	belongs_to_robust_path_final_result.reserve(counter);
	is_robust_path_to_po.reserve(counter);

	for (unsigned int topsorted_iterator = 0; topsorted_iterator <counter; topsorted_iterator++)
	{

			const auto in_degree = boost::in_degree(topsorted_iterator,g);
			if(in_degree==1)// Node must be a NOT-Gate, in this case signal is always robust if there is a transition between the pattern
			{
				auto inEdge = boost::in_edges(topsorted_iterator,g);
				const auto on_input_edge_start_vertex = (*inEdge.first).m_source;
				const auto on_input_edge_end_vertex   = (*inEdge.first).m_target;


				is_robust_or_non_robust_result[on_input_edge_end_vertex].push_back(std::make_pair(on_input_edge_start_vertex,
										circuit_states1[on_input_edge_start_vertex] ^ circuit_states2[on_input_edge_start_vertex])); // in die information fur NOT gates wird doppelt gespeichert
				is_robust_result[on_input_edge_end_vertex].push_back(std::make_pair(on_input_edge_start_vertex,
									    circuit_states1[on_input_edge_start_vertex] ^ circuit_states2[on_input_edge_start_vertex]));

			}
			else if(in_degree==2) // 2 Input case
			{
				auto inEdges = boost::in_edges(topsorted_iterator,g);
				const auto on_input_edge_start_vertex = (*inEdges.first).m_source;
				const auto on_input_edge_end_vertex = (*inEdges.first).m_target;        // equals to topsorted_iterator
				++inEdges.first;                                                       // on_input_edge_end_vertex and off_input_edge_end_vertex are the same, because they point to the current vertex
				const auto off_input_edge_start_vertex = (*inEdges.first).m_source;
				const auto off_input_edge_end_vertex = (*inEdges.first).m_target;      // equals to topsorted_iterator

				if(gate_types_vec[off_input_edge_end_vertex] == 1 || gate_types_vec[off_input_edge_end_vertex] == 4)  // AND and NAND
				{
					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
						isRobustOrNonRobustAndNand(circuit_states1[on_input_edge_start_vertex],
												   circuit_states2[on_input_edge_start_vertex],
												   circuit_states2[off_input_edge_start_vertex])));

					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,
						isRobustOrNonRobustAndNand(circuit_states1[off_input_edge_start_vertex],
												   circuit_states2[off_input_edge_start_vertex],
												   circuit_states2[on_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,isRobustAndNand(circuit_states1[on_input_edge_start_vertex],
																										  circuit_states2[on_input_edge_start_vertex],
																										  circuit_states1[off_input_edge_start_vertex],
																										  circuit_states2[off_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,isRobustAndNand(circuit_states1[off_input_edge_start_vertex],
																										  circuit_states2[off_input_edge_start_vertex],
																										  circuit_states1[on_input_edge_start_vertex],
																										  circuit_states2[on_input_edge_start_vertex])));


				}
				else if(gate_types_vec[off_input_edge_end_vertex] == 5 || gate_types_vec[off_input_edge_end_vertex] == 2)   // Or and NOR case
				{
		
					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
						isRobustOrNonRobustOrNor(circuit_states1[on_input_edge_start_vertex],
												   circuit_states2[on_input_edge_start_vertex],
												   circuit_states2[off_input_edge_start_vertex])));
				
					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,
						isRobustOrNonRobustOrNor(circuit_states1[off_input_edge_start_vertex],
												   circuit_states2[off_input_edge_start_vertex],
												   circuit_states2[on_input_edge_start_vertex])));

					
					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,isRobustOrNor(circuit_states1[on_input_edge_start_vertex],
																										  circuit_states2[on_input_edge_start_vertex],
																										  circuit_states1[off_input_edge_start_vertex],
																										  circuit_states2[off_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,isRobustOrNor(circuit_states1[off_input_edge_start_vertex],
																										  circuit_states2[off_input_edge_start_vertex],
																										  circuit_states1[on_input_edge_start_vertex],
																										  circuit_states2[on_input_edge_start_vertex])));

				}
				else 
				{
					std::cout << "Unknown Gate-Type: " << (int)gate_types_vec[off_input_edge_end_vertex] << std::endl;
	
				}
			}
				
	}
	std::cout << "Path Classification Part 1 ended." << std::endl;

	// Iterate in reverse order over the graph (form the POs to the PIs)
	
	
	for (int reverse_top = counter-1; reverse_top >= 0; reverse_top--)
	{
		const auto in_degree = boost::in_degree(reverse_top,transponse_graph);
		if(in_degree==0)   // An output is found
		{
			is_robust_path_to_po[reverse_top] = boost_bitset(bppk).set();
		}
		else //if(in_degree!=0)
		{
			auto result = boost_bitset(bppk);

			const auto inEdges = boost::in_edges(reverse_top,transponse_graph);
			BOOST_FOREACH(const auto i , inEdges)
			{
					const auto on_input_edge_start_vertex0 = i.m_source;
					const auto on_input_edge_end_vertex0   = i.m_target; // equals reverse_top
	
					const auto t = std::find_if(std::begin(is_robust_result[on_input_edge_start_vertex0]),
							                   std::end(is_robust_result[on_input_edge_start_vertex0]),
							                  [&on_input_edge_end_vertex0](const std::pair<unsigned int,boost_bitset> &o) -> bool 
							            { 
											return on_input_edge_end_vertex0 == o.first; 
										});

					if(t!=is_robust_result[on_input_edge_start_vertex0].cend())
						 result |= is_robust_path_to_po[on_input_edge_start_vertex0]  &  (*t).second;

			}
			is_robust_path_to_po[reverse_top] = std::move(result);
			
		}
		
	}
	std::cout << "Path Classification Part 2 ended." << std::endl;

	// The Last Iteration over the circuit, in Topological order
	for (unsigned int topsorted_iterator = 0; topsorted_iterator <counter; topsorted_iterator++)
	{
		const auto in_degree = boost::in_degree(topsorted_iterator,g);
		//std::cout << "Processing Node: " << vertex_names[topsorted_iterator] << " ";
		if(in_degree!=0)
		{
		
			auto temp_final_result = boost_bitset(bppk);

			const auto inEdges1 = boost::in_edges(topsorted_iterator, g);
			BOOST_FOREACH(const auto edge, inEdges1)
			{
				const auto start_node = edge.m_source;
				const auto end_node = edge.m_target;  // is equal to topsorted_iterator
				/*	auto temp = std::make_pair(on_input_edge_start_vertex,is_robust_path_to_po[on_input_edge_start_vertex] & is_robust_path_to_po[topsorted_iterator];*/  

				/* std::cout << "Processing Edge from: " << vertex_names[start_node] << " to "<< "end_node: " << vertex_names[end_node] <<std::endl;*/
				const auto robust_result = std::find_if(std::begin(is_robust_result[topsorted_iterator]),
														std::end(  is_robust_result[topsorted_iterator]),
					[&start_node](const std::pair<unsigned int, boost_bitset> &o) -> bool
				{ return start_node == o.first; });

				auto belongs_to_robust_path_result_start_node = std::find_if(std::begin(belongs_to_robust_path_final_result[start_node]),
																			 std::end(  belongs_to_robust_path_final_result[start_node]),
					[&end_node](const std::pair<unsigned int, boost_bitset> &o) -> bool
				{ return end_node == o.first; });

				if (robust_result == is_robust_result[topsorted_iterator].cend())
					std::cout << "Could not find the value. " << std::endl;

				if (belongs_to_robust_path_result_start_node == belongs_to_robust_path_final_result[start_node].cend())
					belongs_to_robust_path_result_start_node = belongs_to_robust_path_final_result[start_node].begin();


				/* belongs_to_robust_path_final_result[topsorted_iterator] |= temp.second; */
				temp_final_result |= (is_robust_path_to_po[start_node] & (*belongs_to_robust_path_result_start_node).second & (*robust_result).second);

			}
			if(boost::out_degree(topsorted_iterator,g))
			{
				auto const out_edges = boost::out_edges(topsorted_iterator,g);
				//std::cout << "Processing Edges leaving the node: " << std::endl;
				BOOST_FOREACH(const auto edge,out_edges)
				{
					const auto start_node = edge.m_source; // should eual topsorted_iterator
					const auto end_node = edge.m_target;   // end node
			
					const auto t = std::find_if(std::begin(is_robust_result[end_node]),
							                    std::end(is_robust_result[end_node]),
											   [&start_node](const std::pair<unsigned int,boost_bitset> &o) -> bool 
							                   { return start_node == o.first; } );

					if(t==is_robust_result[end_node].cend())
						std::cout << "Could not find the value. " << std::endl;

					belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(end_node,
																						   temp_final_result &
																						   is_robust_path_to_po[end_node] &
																						   (*t).second));
				}
		
			}
			else
			{

			//	std::cout << "No Edges leaving the node. " << std::endl;
				belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(topsorted_iterator, temp_final_result));
			}
		}
		else
		{
			auto o = boost::out_edges(topsorted_iterator, g);
			BOOST_FOREACH(const auto i, o)
			{
				const auto start_node = i.m_source; // is equal to topsorted_iterator
				const auto end_node = i.m_target;  
				belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(end_node, is_robust_path_to_po[topsorted_iterator]));
			}
			
		}
	}
}

// Path Classification together with logic simulation
void path_classification_with_logic(const Graph                           &g,
	const Graph                           &transponse_graph,
	std::vector<boost_bitset> &circuit_states1,
	std::vector<boost_bitset> &circuit_states2,
	const std::vector<char>                   &gate_types_vec,
	const unsigned int                        &counter,
	const boost_bitset::size_type    &bppk,
	const std::vector<std::string>             &vertex_names,
	std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &is_robust_or_non_robust_result,
	std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &is_robust_result,
	std::vector<boost_bitset> &is_robust_path_to_po,
	std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &belongs_to_robust_path_final_result)
{
	std::cout << "Path Classification starts..." << std::endl;
	is_robust_or_non_robust_result.resize(counter);
	is_robust_result.resize(counter);
	belongs_to_robust_path_final_result.resize(counter);
	is_robust_path_to_po.resize(counter);
	boost_bitset result1(bppk);
	boost_bitset result2(bppk);

	for (unsigned int topsorted_iterator = 0; topsorted_iterator <counter; topsorted_iterator++)
	{

		const auto in_degree = boost::in_degree(topsorted_iterator, g);
		auto inEdges = boost::in_edges(topsorted_iterator, g);  // Get an iterator to the incident edges of the current node
		if (in_degree)
		{
			if (in_degree == 1)// Node must be a NOT-Gate, in this case signal is always robust if there is a transition between the pattern
			{

				result1 = ~circuit_states1[(*inEdges.first).m_source];
				result2 = ~circuit_states2[(*inEdges.first).m_source];

				auto inEdge = boost::in_edges(topsorted_iterator, g);
				const auto on_input_edge_start_vertex = (*inEdge.first).m_source;
				const auto on_input_edge_end_vertex = (*inEdge.first).m_target;


				is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
					circuit_states1[on_input_edge_start_vertex] ^ circuit_states2[on_input_edge_start_vertex])); // in die information fur NOT gates wird doppelt gespeichert
				is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
					circuit_states1[on_input_edge_start_vertex] ^ circuit_states2[on_input_edge_start_vertex]));

			}
			else if (in_degree == 2) // 2 Input case
			{
				result1 = circuit_states1[(*inEdges.first).m_source];
				result2 = circuit_states2[(*inEdges.first).m_source];
				++inEdges.first;

				auto inEdges = boost::in_edges(topsorted_iterator, g);
				const auto on_input_edge_start_vertex = (*inEdges.first).m_source;
				const auto on_input_edge_end_vertex = (*inEdges.first).m_target;        // equals to topsorted_iterator
				++inEdges.first;                                                       // on_input_edge_end_vertex and off_input_edge_end_vertex are the same, because they point to the current vertex
				const auto off_input_edge_start_vertex = (*inEdges.first).m_source;
				const auto off_input_edge_end_vertex = (*inEdges.first).m_target;      // equals to topsorted_iterator

				if (gate_types_vec[off_input_edge_end_vertex] == 1 || gate_types_vec[off_input_edge_end_vertex] == 4)  // AND and NAND
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
					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
						isRobustOrNonRobustAndNand(circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex])));

					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,
						isRobustOrNonRobustAndNand(circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex, isRobustAndNand(circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex],
						circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex, isRobustAndNand(circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex],
						circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex])));


				}
				else if (gate_types_vec[off_input_edge_end_vertex] == 5 || gate_types_vec[off_input_edge_end_vertex] == 2)   // Or and NOR case
				{
					for (; inEdges.first != inEdges.second; ++inEdges.first)
					{
						result1 |= circuit_states1[(*inEdges.first).m_source];
						result2 |= circuit_states2[(*inEdges.first).m_source];
					}
					if (gate_types_vec[topsorted_iterator] == 2)
					{
						result1.flip();
						result2.flip();
					}

					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex,
						isRobustOrNonRobustOrNor(circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex])));

					is_robust_or_non_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex,
						isRobustOrNonRobustOrNor(circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex])));


					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(on_input_edge_start_vertex, isRobustOrNor(circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex],
						circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex])));

					is_robust_result[on_input_edge_end_vertex].emplace_back(std::make_pair(off_input_edge_start_vertex, isRobustOrNor(circuit_states1[off_input_edge_start_vertex],
						circuit_states2[off_input_edge_start_vertex],
						circuit_states1[on_input_edge_start_vertex],
						circuit_states2[on_input_edge_start_vertex])));

				}
				else if (gate_types_vec[topsorted_iterator] == 3)//if(*str_iter=="XOR")
				{
					for (; inEdges.first != inEdges.second; ++inEdges.first)
					{
						result1 ^= circuit_states1[(*inEdges.first).m_source];
						result2 ^= circuit_states2[(*inEdges.first).m_source];
					}

					std::cout << "Unknown Gate-Type, XOR: " << (int)gate_types_vec[off_input_edge_end_vertex] << std::endl;

				}
			}
			circuit_states1[topsorted_iterator] = std::move(result1);
			circuit_states2[topsorted_iterator] = std::move(result2);
		}
	}
	std::cout << "Path Classification Part 1 ended." << std::endl;

	// Iterate in reverse order over the graph (form the POs to the PIs)


	for (int reverse_top = counter - 1; reverse_top >= 0; reverse_top--)
	{
		const auto in_degree = boost::in_degree(reverse_top, transponse_graph);
		if (in_degree == 0)   // An output is found
		{
			is_robust_path_to_po[reverse_top] = boost_bitset(bppk).set();
		}
		else //if(in_degree!=0)
		{
			auto result = boost_bitset(bppk);

			const auto inEdges = boost::in_edges(reverse_top, transponse_graph);
			BOOST_FOREACH(const auto i, inEdges)
			{
				const auto on_input_edge_start_vertex0 = i.m_source;
				const auto on_input_edge_end_vertex0 = i.m_target; // equals reverse_top

				const auto t = std::find_if(std::cbegin(is_robust_result[on_input_edge_start_vertex0]),
					std::cend(is_robust_result[on_input_edge_start_vertex0]),
					[&on_input_edge_end_vertex0](const std::pair<unsigned int, boost_bitset> &o) -> bool
				{
					return on_input_edge_end_vertex0 == o.first;
				});

				if (t != is_robust_result[on_input_edge_start_vertex0].cend())
					result |= is_robust_path_to_po[on_input_edge_start_vertex0] & (*t).second;

			}
			is_robust_path_to_po[reverse_top] = std::move(result);

		}

	}
	std::cout << "Path Classification Part 2 ended." << std::endl;

	// The Last Iteration over the circuit, in Topological order
	for (unsigned int topsorted_iterator = 0; topsorted_iterator <counter; topsorted_iterator++)
	{
		const auto in_degree = boost::in_degree(topsorted_iterator, g);
		//std::cout << "Processing Node: " << vertex_names[topsorted_iterator] << " ";
		if (in_degree != 0)
		{

			auto temp_final_result = boost_bitset(bppk);

			const auto inEdges1 = boost::in_edges(topsorted_iterator, g);
			BOOST_FOREACH(const auto edge, inEdges1)
			{
				const auto start_node = edge.m_source;
				const auto end_node = edge.m_target;  // is equal to topsorted_iterator
				/*	auto temp = std::make_pair(on_input_edge_start_vertex,is_robust_path_to_po[on_input_edge_start_vertex] & is_robust_path_to_po[topsorted_iterator];*/

				/* std::cout << "Processing Edge from: " << vertex_names[start_node] << " to "<< "end_node: " << vertex_names[end_node] <<std::endl;*/
				const auto robust_result = std::find_if(std::cbegin(is_robust_result[topsorted_iterator]),
					std::cend(is_robust_result[topsorted_iterator]),
					[&start_node](const std::pair<unsigned int, boost_bitset> &o) -> bool
				{ return start_node == o.first; });

				auto belongs_to_robust_path_result_start_node = std::find_if(std::cbegin(belongs_to_robust_path_final_result[start_node]),
					std::cend(belongs_to_robust_path_final_result[start_node]),
					[&end_node](const std::pair<unsigned int, boost_bitset> &o) -> bool
				{ return end_node == o.first; });

				if (robust_result == is_robust_result[topsorted_iterator].cend())
					std::cout << "Could not find the value. " << std::endl;

				if (belongs_to_robust_path_result_start_node == belongs_to_robust_path_final_result[start_node].cend())
					belongs_to_robust_path_result_start_node = belongs_to_robust_path_final_result[start_node].begin();


				/* belongs_to_robust_path_final_result[topsorted_iterator] |= temp.second; */
				temp_final_result |= (is_robust_path_to_po[start_node] & (*belongs_to_robust_path_result_start_node).second & (*robust_result).second);

			}
			if (boost::out_degree(topsorted_iterator, g))
			{
				auto const out_edges = boost::out_edges(topsorted_iterator, g);
				//std::cout << "Processing Edges leaving the node: " << std::endl;
				BOOST_FOREACH(const auto edge, out_edges)
				{
					const auto start_node = edge.m_source; // should eual topsorted_iterator
					const auto end_node = edge.m_target;   // end node

					const auto t = std::find_if(std::cbegin(is_robust_result[end_node]),
						std::cend(is_robust_result[end_node]),
						[&start_node](const std::pair<unsigned int, boost_bitset> &o) -> bool
					{ return start_node == o.first; });

					if (t == is_robust_result[end_node].cend())
						std::cout << "Could not find the value. " << std::endl;

					belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(end_node,
						temp_final_result &
						is_robust_path_to_po[end_node] &
						(*t).second));
				}

			}
			else
			{

				//	std::cout << "No Edges leaving the node. " << std::endl;
				belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(topsorted_iterator, temp_final_result));
			}
		}
		else
		{
			auto o = boost::out_edges(topsorted_iterator, g);
			BOOST_FOREACH(const auto i, o)
			{
				const auto start_node = i.m_source; // is equal to topsorted_iterator
				const auto end_node = i.m_target;
				belongs_to_robust_path_final_result[topsorted_iterator].emplace_back(std::make_pair(end_node, is_robust_path_to_po[topsorted_iterator]));
			}

		}
	}
}

void logicsimulation(unsigned int length,
	std::vector<unsigned int> top_order,
	std::vector<unsigned int> reverse_top_order,
	std::vector<unsigned int> gate_type,
	std::vector<std::pair<unsigned int,unsigned int>> edges,
	std::vector<std::pair<unsigned int, unsigned int>> out_edges,
	std::vector<unsigned int> &p1,
	std::vector<unsigned int> &p2,
	std::vector<std::pair<unsigned int,unsigned int>> &robust_result,
	std::vector<std::pair<unsigned int, unsigned int >> &is_robust_or_non_robust_result,
	std::vector<unsigned int> &is_robust_path_to_po,
	std::vector<std::pair<unsigned int, unsigned int>> &belongs_to_robust_path_final_result,
	const std::vector<std::string>             &vertex_names)
{
	unsigned int output_offset = length*0;
	//printf("%i ",tid);
	int size = length;
	unsigned int  result1 = 0;
	unsigned int  result2 = 0;
	//int pattern_counter = 0;
	unsigned int pattern1_1;
	unsigned int pattern1_2;
	unsigned int pattern2_1;
	unsigned int pattern2_2;
	unsigned int isRobust_1;
	unsigned int isRobust_2;
	unsigned int isRobustOrNonRobust_l;
	unsigned int isRobustOrNonRobust_2;
	for (int i = 0; i<size; i++)
	{
		unsigned int c_node = top_order[i];
		result1 = 0;
		result2 = 0;
		unsigned int g = gate_type[c_node];
		if (g)   // All other gates, if node is an input the else is executed, the gate type can also be used to determine the in_degree
		{

			pattern1_1 = p1[output_offset + edges[c_node].first];
			pattern1_2 = p1[output_offset + edges[c_node].second];
			pattern2_1 = p2[output_offset + edges[c_node].first];
			pattern2_2 = p2[output_offset + edges[c_node].second];

			if (g == 1)
			{

				result1 = ~(pattern1_1 & pattern1_2);
				result2 = ~(pattern2_1 & pattern2_2);
			}
			else if (g == 2)
			{

				result1 = ~(pattern1_1 | pattern1_2);
				result2 = ~(pattern2_1 | pattern2_2);
			}
			else if (g == 4)
			{

				result1 = pattern1_1 & pattern1_2;
				result2 = pattern2_1 & pattern2_2;
			}
			else if (g == 5)
			{

				result1 = pattern1_1 | pattern1_2;
				result2 = pattern2_1 | pattern2_2;
			}
			else if (g == 6)
			{
				result1 = ~pattern1_1;
				result2 = ~pattern2_1;
			}

			if (g == 1 || g == 4) // AND NAND case
			{
				 // isRobustAND NAND
				isRobust_1 = ((~pattern1_1) & pattern1_2 & pattern2_2) | (pattern1_1& (~pattern1_2) & pattern2_1 & pattern2_2);
				isRobust_2 = ((~pattern2_1) & pattern2_2 & pattern1_2) | (pattern2_1& (~pattern2_2) & pattern1_1 & pattern1_2);

				isRobustOrNonRobust_l = (pattern1_1 ^pattern1_2) & pattern2_2;
				isRobustOrNonRobust_2 = (pattern2_1 ^pattern2_2) & pattern1_2;
			}
			else if (g == 5 || g == 2) // Or Nor
			{
				isRobust_1 = (pattern1_1 & (~pattern1_2) & (~pattern2_2)) | ((~pattern1_1) & pattern1_2 & (~pattern2_1) & (~pattern2_2));
				isRobust_2 = (pattern2_1 & (~pattern2_2) & (~pattern1_2)) | ((~pattern2_1) & pattern2_2 & (~pattern1_1) & (~pattern1_2));


				isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & (~pattern2_2);
				isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & (~pattern1_2);
			}
			else if (g == 6)
			{
				isRobust_1 = pattern1_1 ^ pattern2_1;
				isRobust_2 = pattern1_1 ^ pattern2_1;

				isRobustOrNonRobust_l = pattern1_1 ^ pattern2_1;
				isRobustOrNonRobust_2 = pattern1_1 ^ pattern2_1;
			}
			
			p1[output_offset + c_node] = result1;
			p2[output_offset + c_node] = result2;

			robust_result[output_offset + c_node] = std::make_pair(isRobust_1, isRobust_2);
			is_robust_or_non_robust_result[output_offset + c_node] = std::make_pair(isRobustOrNonRobust_l, isRobustOrNonRobust_2);

		}
		else
		{
			p1[output_offset + c_node] = boost_bitset(generate4bit()).to_ulong();
			p2[output_offset + c_node] = boost_bitset(generate4bit()).to_ulong();
		}
	}

	/* traverse over the graph in reverse topological order */
	for (int i = 0; i < size; i++)
	{
		unsigned int c_node = reverse_top_order[i];
		if ((out_edges[c_node].first != 0) && (out_edges[c_node].second != 0)) // Node is not an output
		{
			unsigned int result = 0;
			
			unsigned int r1 = 0;
			unsigned int path_po = 0;
			if (edges[out_edges[c_node].first].first == c_node)
			{
				r1      = robust_result[output_offset + out_edges[c_node].first].first;
				path_po = is_robust_path_to_po[output_offset + out_edges[c_node].first];
			}
			else
			{
				r1      = robust_result[output_offset + out_edges[c_node].first].second;
				path_po = is_robust_path_to_po[output_offset + out_edges[c_node].first];
			}

			result |= (r1 & path_po);

			if (edges[out_edges[c_node].second].first == c_node)
			{
				r1 =       robust_result[output_offset + out_edges[c_node].second].first;
				path_po = is_robust_path_to_po[output_offset + out_edges[c_node].second];
			}
			else
			{
				r1 =     robust_result[output_offset + out_edges[c_node].second].second;
				path_po = is_robust_path_to_po[output_offset + out_edges[c_node].second];
			}
			result |= (r1 & path_po);
			is_robust_path_to_po[output_offset + c_node] = result;
		}
		else
		{
			is_robust_path_to_po[output_offset + c_node] = ~0; // set the bits of all outputs to 1
		}// Node is an output because no out edges
		//case of output set all bits to "1"
	}

	/* Last Interation over the circuit in topological order*/
	for (int i = 0; i < size; i++)
	{
		unsigned int c_node = top_order[i];
		unsigned int g = gate_type[c_node];
		if (g)
		{
			/* Check all input edges */
			unsigned int temp_final_result = 0;
			unsigned int btrofr = 0;
			/* robust_result[output_offset + c_node].first; */

			if (out_edges[edges[c_node].first].first == c_node)
				btrofr = belongs_to_robust_path_final_result[output_offset + edges[c_node].first].first;
			else
				btrofr = belongs_to_robust_path_final_result[output_offset + edges[c_node].first].second;

			temp_final_result |= is_robust_path_to_po[output_offset + edges[c_node].first] & robust_result[output_offset + c_node].first & btrofr;

			/* robust_result[output_offset + c_node].second; */
			if (out_edges[edges[c_node].second].first == c_node)
				btrofr = belongs_to_robust_path_final_result[output_offset + edges[c_node].second].first;
			else
				btrofr = belongs_to_robust_path_final_result[output_offset + edges[c_node].second].second;

			temp_final_result |= is_robust_path_to_po[output_offset + edges[c_node].second] & robust_result[output_offset + c_node].second & btrofr;

			/* Check all output edges */
			if ((out_edges[c_node].first != 0) && (out_edges[c_node].second != 0))
			{
				unsigned int result1 = 0;
				unsigned int result2 = 0;
				
				unsigned int t = 0;
				/*out_edges[c_node].first; */
				if (edges[out_edges[c_node].first].first == c_node)
					t = robust_result[output_offset + out_edges[c_node].first].first;
				else
					t = robust_result[output_offset + out_edges[c_node].first].second;
				result1 = t & temp_final_result & is_robust_path_to_po[output_offset + out_edges[c_node].first];

				if (edges[out_edges[c_node].second].first == c_node)
					t = robust_result[output_offset + out_edges[c_node].second].first;
				else
					t = robust_result[output_offset + out_edges[c_node].second].second;
				result2 = t & temp_final_result & is_robust_path_to_po[output_offset + out_edges[c_node].second];

				belongs_to_robust_path_final_result[output_offset + c_node] = std::make_pair(result1,result2);

			}
			else
			{
				belongs_to_robust_path_final_result[output_offset + c_node] = 
					std::make_pair(temp_final_result, temp_final_result);
			}
		}
		else
		{
			belongs_to_robust_path_final_result[output_offset + c_node].first = is_robust_path_to_po[output_offset + c_node];
			belongs_to_robust_path_final_result[output_offset + c_node].second = is_robust_path_to_po[output_offset + c_node];
		}
	}
}


void logicsimulationv2(
	const unsigned int number_of_threads,
	const unsigned int length,
	const std::vector<cl_uint> top_order,
	const std::vector<cl_uint> reverse_top_order,
	const std::vector<cl_uint> gate_type,
	const std::vector<std::pair<cl_uint, cl_uint>> edges,
	const std::vector<std::pair<cl_uint, cl_uint>> out_edges,
	std::vector<cl_uint> &p1,
	std::vector<cl_uint> &p2,
	std::vector<cl_uint> &robust_result1,
	std::vector<cl_uint> &robust_result2,
	std::vector<cl_uint> &is_robust_or_non_robust_result1,
	std::vector<cl_uint> &is_robust_or_non_robust_result2,
	std::vector<cl_uint> &is_robust_path_to_po,
	std::vector<cl_uint> &belongs_to_robust_path_final_result1,
	std::vector<cl_uint> &belongs_to_robust_path_final_result2)
{
	for (std::size_t c = 0; c < number_of_threads; c++)
	{
		auto output_offset = length * c;
	
		int size = length;
		unsigned int  result1 = 0;
		unsigned int  result2 = 0;
		//int pattern_counter = 0;
		unsigned int pattern1_1;
		unsigned int pattern1_2;
		unsigned int pattern2_1;
		unsigned int pattern2_2;
		unsigned int isRobust_1;
		unsigned int isRobust_2;
		unsigned int isRobustOrNonRobust_l;
		unsigned int isRobustOrNonRobust_2;
		for (int i = 0; i < size; i++)
		{
			unsigned int c_node = top_order[i];
			result1 = 0;
			result2 = 0;
			unsigned int g = gate_type[c_node];
			if (g)   // All other gates, if node is an input the else is executed, the gate type can also be used to determine the in_degree
			{

				pattern1_1 = p1[output_offset + edges[c_node].first];
				pattern1_2 = p1[output_offset + edges[c_node].second];
				pattern2_1 = p2[output_offset + edges[c_node].first];
				pattern2_2 = p2[output_offset + edges[c_node].second];

				if (g == 1)
				{

					result1 = ~(pattern1_1 & pattern1_2);
					result2 = ~(pattern2_1 & pattern2_2);
				}
				else if (g == 2)
				{

					result1 = ~(pattern1_1 | pattern1_2);
					result2 = ~(pattern2_1 | pattern2_2);
				}
				else if (g == 4)
				{

					result1 = pattern1_1 & pattern1_2;
					result2 = pattern2_1 & pattern2_2;
				}
				else if (g == 5)
				{

					result1 = pattern1_1 | pattern1_2;
					result2 = pattern2_1 | pattern2_2;
				}
				else if (g == 6)
				{
					result1 = ~pattern1_1;
					result2 = ~pattern2_1;
				}

				if (g == 1 || g == 4) // AND NAND case
				{
					// isRobustAND NAND
					isRobust_1 = ((~pattern1_1) & pattern1_2 & pattern2_2) | (pattern1_1& (~pattern1_2) & pattern2_1 & pattern2_2);
					isRobust_2 = ((~pattern2_1) & pattern2_2 & pattern1_2) | (pattern2_1& (~pattern2_2) & pattern1_1 & pattern1_2);

					isRobustOrNonRobust_l = (pattern1_1 ^pattern1_2) & pattern2_2;
					isRobustOrNonRobust_2 = (pattern2_1 ^pattern2_2) & pattern1_2;
				}
				else if (g == 5 || g == 2) // Or Nor
				{
					isRobust_1 = (pattern1_1 & (~pattern1_2) & (~pattern2_2)) | ((~pattern1_1) & pattern1_2 & (~pattern2_1) & (~pattern2_2));
					isRobust_2 = (pattern2_1 & (~pattern2_2) & (~pattern1_2)) | ((~pattern2_1) & pattern2_2 & (~pattern1_1) & (~pattern1_2));


					isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & (~pattern2_2);
					isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & (~pattern1_2);
				}
				else if (g == 6)
				{
					isRobust_1 = pattern1_1 ^ pattern2_1;
					isRobust_2 = pattern1_1 ^ pattern2_1;

					isRobustOrNonRobust_l = pattern1_1 ^ pattern2_1;
					isRobustOrNonRobust_2 = pattern1_1 ^ pattern2_1;
				}

				p1[output_offset + c_node] = result1;
				p2[output_offset + c_node] = result2;

				robust_result1[output_offset + c_node] = isRobust_1;
				robust_result2[output_offset + c_node] = isRobust_2;
				is_robust_or_non_robust_result1[output_offset + c_node] = isRobustOrNonRobust_l;
				is_robust_or_non_robust_result2[output_offset + c_node] = isRobustOrNonRobust_2;

			}
			else
			{
				p1[output_offset + c_node] = boost_bitset(generate4bit()).to_ulong();
				p2[output_offset + c_node] = boost_bitset(generate4bit()).to_ulong();
			}
		}

		/* traverse over the graph in reverse topological order */
		for (int i = 0; i < size; i++)
		{
			unsigned int c_node = reverse_top_order[i];
			if ((out_edges[c_node].first != 0) && (out_edges[c_node].second != 0)) // Node is not an output
			{
				unsigned int result = 0;

				unsigned int r1 = 0;
				unsigned int path_po = 0;
				if (edges[out_edges[c_node].first].first == c_node)
				{
					r1 = robust_result1[output_offset + out_edges[c_node].first];
					path_po = is_robust_path_to_po[output_offset + out_edges[c_node].first];
				}
				else
				{
					r1 = robust_result2[output_offset + out_edges[c_node].first];
					path_po = is_robust_path_to_po[output_offset + out_edges[c_node].first];
				}

				result |= (r1 & path_po);

				if (edges[out_edges[c_node].second].first == c_node)
				{
					r1 = robust_result1[output_offset + out_edges[c_node].second];
					path_po = is_robust_path_to_po[output_offset + out_edges[c_node].second];
				}
				else
				{
					r1 = robust_result2[output_offset + out_edges[c_node].second];
					path_po = is_robust_path_to_po[output_offset + out_edges[c_node].second];
				}
				result |= (r1 & path_po);
				is_robust_path_to_po[output_offset + c_node] = result;
			}
			else
			{
				is_robust_path_to_po[output_offset + c_node] = ~0; // set the bits of all outputs to 1
			}// Node is an output because no out edges
			//case of output set all bits to "1"
		}

		/* Last Interation over the circuit in topological order*/
		for (int i = 0; i < size; i++)
		{
			unsigned int c_node = top_order[i];
			unsigned int g = gate_type[c_node];
			if (g)
			{
				/* Check all input edges */
				unsigned int temp_final_result = 0;
				unsigned int btrofr = 0;
				/* robust_result[output_offset + c_node].first; */

				if (out_edges[edges[c_node].first].first == c_node)
					btrofr = belongs_to_robust_path_final_result1[output_offset + edges[c_node].first];
				else
					btrofr = belongs_to_robust_path_final_result2[output_offset + edges[c_node].first];

				temp_final_result |= is_robust_path_to_po[output_offset + edges[c_node].first] & robust_result1[output_offset + c_node] & btrofr;

				/* robust_result[output_offset + c_node].second; */
				if (out_edges[edges[c_node].second].first == c_node)
					btrofr = belongs_to_robust_path_final_result1[output_offset + edges[c_node].second];
				else
					btrofr = belongs_to_robust_path_final_result2[output_offset + edges[c_node].second];

				temp_final_result |= is_robust_path_to_po[output_offset + edges[c_node].second] & robust_result2[output_offset + c_node] & btrofr;

				/* Check all output edges */
				if ((out_edges[c_node].first != 0) && (out_edges[c_node].second != 0))
				{
					unsigned int result1 = 0;
					unsigned int result2 = 0;

					unsigned int t = 0;
					/*out_edges[c_node].first; */
					if (edges[out_edges[c_node].first].first == c_node)
						t = robust_result1[output_offset + out_edges[c_node].first];
					else
						t = robust_result2[output_offset + out_edges[c_node].first];
					result1 = t & temp_final_result & is_robust_path_to_po[output_offset + out_edges[c_node].first];

					if (edges[out_edges[c_node].second].first == c_node)
						t = robust_result1[output_offset + out_edges[c_node].second];
					else
						t = robust_result2[output_offset + out_edges[c_node].second];
					result2 = t & temp_final_result & is_robust_path_to_po[output_offset + out_edges[c_node].second];

					belongs_to_robust_path_final_result1[output_offset + c_node] = result1;
					belongs_to_robust_path_final_result2[output_offset + c_node] = result2;

				}
				else
				{
					belongs_to_robust_path_final_result1[output_offset + c_node] = temp_final_result;
					belongs_to_robust_path_final_result2[output_offset + c_node] = temp_final_result;
				}
			}
			else
			{
				belongs_to_robust_path_final_result1[output_offset + c_node] = is_robust_path_to_po[output_offset + c_node];
				belongs_to_robust_path_final_result2[output_offset + c_node] = is_robust_path_to_po[output_offset + c_node];
			}
		}
	}
}
#endif