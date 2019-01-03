/*
   This file contains the create_graph function. With the help of this funtions the graph can be created from a 
   previously parsed netlist which is in the "circuit &c" parameter.
   Even though the function is very big the algorithm is very easy. It iterates over all gates and created an edge
   between every input and output(s). In addition to creating the graph also some meta information is collected for
   debugging 
*/


#ifndef GRAPH_CREATOR
#define GRAPH_CREATOR

#include <boost/algorithm/string.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <unordered_map>

#include "file_parsers.h"
#include "hash_maps.h"
#include "helper_functions.h"
#include "graph_structure.h"

Graph create_graph(const circuit &c,
	std::vector<std::string> &node_names,
	std::vector<std::string> &gate_type_str,
	std::vector<char> &gate_types_char,
	std::vector<std::pair<bool, unsigned int>> &constant_input_or_wires,
	unsigned int &final_number_of_nodes,
	my_map &io_map,
	dff_map &ff_info_map)
{
	const std::size_t estimated_size = c.circuit_gates.size() + c.circuit_gates.size()/4; 
	std::vector<char>gate_types_vec;                                                      // the gate type as number e.g 0 for NAND
	std::vector<std::string>vertex_names;                                                 // the name of the input or output
	std::vector<std::string>gate_type;                                                    // the gate type as string e.g. AND2X, same like gate_types_vec
	std::vector<std::pair<unsigned int, unsigned int>>edges;                              // this vector holds the edges

	gate_types_vec.reserve(estimated_size);                                               // in the following lines some memory is reserved, even though
	vertex_names.reserve(estimated_size);                                                 // the exact amount is unknown. If the size is too small C++
	gate_type.reserve(estimated_size);                                                    // will automaticall increase the size
	edges.reserve(2 * estimated_size); 

	unsigned int counter = 0;

	const std::vector<std::string> outputs1 = { "Q", "QN", "ZN", "Z" };                  // These are the outputs the algorithm is looking for
	for (const auto &gate_iterator : c.circuit_gates)                                    // the loop is iterating over all gates in the netlist
	{
		if (gate_iterator.gate_type.find("DFF") == std::string::npos)                // Gate is not a FlipFlop, FlipFlops are skipped
		{
			std::vector<std::pair<std::string,std::string>>found_outputs;            
			for (const auto &i : gate_iterator.ipc)                              // this loop finds the outputs of each gate and saves them in "found_outputs"
			{
				          
				const auto found_out = std::find(std::cbegin(outputs1), std::cend(outputs1), boost::to_upper_copy(i.first));

				if (found_out != std::cend(outputs1))
					found_outputs.emplace_back(i.second,i.first);
			}
		
			for (const auto &i : gate_iterator.ipc)    // in this loop edges between the inputs and outputs are created
			{
				if (std::find(std::cbegin(outputs1), std::cend(outputs1), i.first) == std::cend(outputs1)) // if the pin it not an output process it
				{
					std::string str = (i.second.find("1\'") != std::string::npos || i.second.find("0\'") != std::string::npos) ? 
									  std::string(i.second.cbegin() + 2, i.second.cend()) : 
									  i.second; // because inputs can be constant e.g. 0' we must check first and cut this part if it true

			
					const auto r = io_map.emplace(str, std::make_pair(i.first, counter));   // insert the current signal into a hashmap
					if (r.second)                                                           // if the node was not present in the hash map it was created and in the if-body some meta information about the node gate-type etc. will be saved, if the node already exists the body of this if-block will not be executed
					{
						if (i.second.find("1\'") != std::string::npos || i.second.find("0\'") != std::string::npos)
						{
							constant_input_or_wires.emplace_back(i.second[0]=='1' ? true : false, counter);
							vertex_names.emplace_back(std::move(str));
						}
						else vertex_names.emplace_back(str); 
						gate_types_vec.emplace_back(gateTypeToInteger((gate_iterator).gate_type));
						gate_type.emplace_back(gate_iterator.gate_type);
						counter++; // that means the hash-entry did noty exist and it was created
					}

					// now create an edge between the current input and the output(s)
					for (const auto &o : found_outputs)  // found_outputs contains the previoulsy found outputs of the current gate
					{
						const auto output_node = io_map.emplace(o.first, std::make_pair(o.second, counter)); // check if the current output already has an entry in the hash table, if no the if-body will be executed in order to record the meta information
						if (output_node.second)                                                                               // if the output node did not exist
						{

							vertex_names.emplace_back(o.first);
							if (o.second.find("1\'") != std::string::npos)constant_input_or_wires.emplace_back(true, counter);
							else if (o.second.find("0\'") != std::string::npos)constant_input_or_wires.emplace_back(false, counter);
							gate_types_vec.emplace_back(gateTypeToInteger((gate_iterator).gate_type));
							gate_type.emplace_back(gate_iterator.gate_type);
							counter++; // that means the hash-entry did noty exist and it was created
						}
						// if the entry already exists in the hash-table and edge is created
						gate_types_vec[(*output_node.first).second.second] = gateTypeToInteger((gate_iterator).gate_type);
						gate_type[(*output_node.first).second.second] = gate_iterator.gate_type;
						edges.emplace_back((*r.first).second.second, (*output_node.first).second.second);
								
					}
				}
			}
		}
	}

    // print some information about the circuit
	std::cout << "Gates in circuit: " << c.circuit_gates.size() << std::endl;
	std::cout << "Size of Vertex names: " << vertex_names.size() << std::endl;


	std::cout << "Number of Vertices: " << counter << " Number of Edges: " << edges.size() << std::endl;

	Graph g(edges.begin(), edges.end(), counter);   // here the graph is created by adding the edges to it

	// if the next three lines are uncommented the graph is checked for cycles
	/*bool has_cycle = false;
	cycle_detector vis(has_cycle);
	boost::depth_first_search(g, visitor(vis));
	std::cout << "The graph has a cycle? " << (has_cycle ? "Yes" : "No") << std::endl;*/

	// the next four lines create the topological order by using the boost::topological_sort function from the graph library
	// because this function returns the toplogical order but in reverse order (not reverse toplogical order), that means the
	// inputs are found at the end of the array, for this reason std::reverse is called to reverse the array
	std::vector<unsigned int> top_sorted_graph;
	top_sorted_graph.reserve(counter);
	boost::topological_sort(g, std::back_inserter(top_sorted_graph)); 
	std::reverse(std::begin(top_sorted_graph), std::end(top_sorted_graph));

	// Create a levelized toplogical order of the circuit, if the next three lines are commented the normal toplogical order will be used, otherwise a levelized toplogical order
	std::vector<unsigned int> lvl_ordered;
	levelice_circuit(g, top_sorted_graph, counter, lvl_ordered);
	std::copy(std::cbegin(lvl_ordered), std::cend(lvl_ordered), std::begin(top_sorted_graph));
	std::cout << "Topological Sort finished." << std::endl;

	// From here until the end of this function the meta data for each input/output is ordered according to the toplogical order
	// This loop is iterating over the maps each node in the toplogical order to its array position i
	std::unordered_map<unsigned int, unsigned int>top_order_map;
	for (unsigned int i = 0; i < counter; i++)
	{
		const auto j = top_sorted_graph[i];
		top_order_map.insert(std::make_pair(j, i));
	}

	node_names.reserve(counter);
	gate_type_str.reserve(counter);
	gate_types_char.reserve(counter);

	// This loop orders the meta information of each node topological
	for (unsigned int i = 0; i < counter; i++)
	{
		const auto j = top_sorted_graph[i];
		gate_types_char.emplace_back(gate_types_vec[j]);
		node_names.emplace_back(vertex_names[j]);
		gate_type_str.emplace_back(gate_type[j]);
	}

	// This loop iterates over the "old" edges and creates "new" edges which use
	// the topological order
	for (auto& e : edges)
	{
		auto i1 = top_order_map.find(e.first);
		auto i2 = top_order_map.find(e.second);
		if (i1 == std::cend(top_order_map) || i2 == std::cend(top_order_map))
			std::cout << "Could not find edge in Hash-Table" << std::endl;
		else
		{
			e.first = (*i1).second;
			e.second = (*i2).second;
		}
	}
	// This loop updates the constant wires information e.g 0' or 1'
	for (auto &c : constant_input_or_wires)
	{
		auto i = top_order_map.find(c.second);
		c.second = (*i).second;
	}
	// this loop updates information in a hash table which is used when using the stilefile
	for (auto &i : io_map)
	{
		auto i1 = top_order_map.find(i.second.second);
		if (i1 == std::cend(top_order_map))
			std::cout << "Could not find edge in Hash-Table" << std::endl;
		else
			i.second.second = (*i1).second;
	}
	final_number_of_nodes = counter;
	g = Graph(edges.begin(), edges.end(), counter); // The final graph is created and returned
	return g;
}

#endif
