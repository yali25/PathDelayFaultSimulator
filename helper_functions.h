/*
  This file defines several functions which may help to develop
  VLSI algrithms quickly
*/

#ifndef HELPER_FUNCTIONS
#define HELPER_FUNCTIONS
#include <vector>
#include <algorithm>
#include <string>
#include <random>
#include <fstream>
#include <boost/graph/transpose_graph.hpp>
#include <boost/dynamic_bitset.hpp>
#include <CL/cl.hpp>
#include "graph_structure.h"
#include "hash_maps.h"

// with the help of this function the levelized topological ordered graph
// can be obtained, this first parameter is the graph, followed by the topological order,
// counter represents the total number of nodes in the graph
// at the end of the function lvl_ordered will contain the levelized order of the graph
void levelice_circuit(const Graph &g,
					  const std::vector<unsigned int> &top_sorted_graph,
					  const unsigned &counter,
					  std::vector<unsigned int> &lvl_ordered)
{
	std::vector<unsigned int> lvl_info; 
	std::vector<std::vector<unsigned int>> levelized_circuit;
	lvl_info.resize(counter);
	levelized_circuit.resize(1);
	levelized_circuit.push_back(std::vector<unsigned int>());
	for (const auto topsorted_iterator : top_sorted_graph)     // this loop is iterating over the topological ordered graph
	{
			if(boost::in_degree(topsorted_iterator,g)==0)     // if the in_degree is 0, it means an input is found, in this case put the node in the first level (0)
			{
				lvl_info[topsorted_iterator] = 0;                    // for the current node we save the level because later we must use it again for other nodes which have a connection to the current node          
				levelized_circuit[0].push_back(topsorted_iterator);  // 0 is the first level
			}
			else 
			{   // if the in_degree != 0 it means the current node is not an input 
				auto inEdges = boost::in_edges(topsorted_iterator,g);
				std::vector<unsigned int>inEdgesLvl;				
				for(; inEdges.first != inEdges.second; ++inEdges.first)    // this loop obtains all levels of the parent nodes
				{
						inEdgesLvl.push_back(lvl_info[(*inEdges.first).m_source]);  // save the levels of the parent node 
				}
				const unsigned int level =  *std::max_element(inEdgesLvl.cbegin(),inEdgesLvl.cend())+1; // now find the maximum level of the parent node and add 1 to it for the current node
				lvl_info[topsorted_iterator] = level;                                                   // assign the calculated level to the current node
				if(levelized_circuit.size()-1 < level)                                                  // it may happen that the levelized_circuit data structure does not contain the current level, in this case add a new level
				{
					std::vector<unsigned int>vec;                                                       // the new level
                    vec.push_back(topsorted_iterator);                                                  // the current node is on that level
					levelized_circuit.push_back(vec);                                                    
				}
				else
					levelized_circuit[lvl_info[topsorted_iterator]].push_back(topsorted_iterator);     // level already exsists, no need to create a new level
				lvl_info[topsorted_iterator] = level;

			}
		}
		// at the end check if all nodes have a level, and copy from the 2D data-structure to a 1D data-structure 
		std::cout << "Verify that all Nodes have a lvl: "; 
		std::size_t all_lvl_sizes = 0;
		unsigned int lvl = 1;
		std::ofstream f("leon3mp_lvl_info.txt");
		for(const auto lvl_iter :  levelized_circuit)
		{
			f << lvl << " " << (lvl_iter).size() << std::endl;
			lvl++;
			all_lvl_sizes += (lvl_iter).size();
			std::copy(std::cbegin(lvl_iter), std::cend(lvl_iter), std::back_inserter(lvl_ordered));  // copy each level to the 1D data structure
		}
		std::cout << all_lvl_sizes << "==" << counter << std::endl;
}


// This class is used by the Boost graph library to find cycles inside
// the graph, the program uses this class only once, even though it is
// not necessary because the graph can not have any cycle because of the
// way the graph is created avoids cycles.
 struct cycle_detector : public boost::dfs_visitor<>
 {
   cycle_detector( bool& has_cycle): _has_cycle(has_cycle) { }

    template <class Edge, class Graph>
    void back_edge(Edge, Graph&) {
      _has_cycle = true;
    }
  protected:
    bool& _has_cycle;
  };


void stilfile_check(const stilfile_datastruct &s, 
				  const my_map &io_map, 
				  const dff_map &ff_info_map,
				  std::vector<boost_bitset > &circuit_states1,
				  bitset_map &po_expected_values)
{
	std::vector<scan_chain>::const_iterator scan_chain_it;
	std::vector<std::vector<std::size_t>>::iterator k;
	std::vector<std::vector<std::size_t>>scan_chain_struc_lookup;
	scan_chain_struc_lookup.resize(s.scan_chains.size());

	for(scan_chain_it = s.scan_chains.cbegin(),
	    k = scan_chain_struc_lookup.begin();
	   scan_chain_it!=s.scan_chains.cend();scan_chain_it++,k++)
    {
	   for(auto j=(*scan_chain_it).scan_structure.rbegin();j!=(*scan_chain_it).scan_structure.rend();j++)
	   {
		   dff_map::const_iterator lr = ff_info_map.find((*j).ff_name);
		   if(lr == ff_info_map.end())
		   {
			   std::cout << "Could not find element: " << (*j).ff_name << " in hashtable." << std::endl;
		   }
		   else
		   {
			   my_map::const_iterator m = io_map.find((*lr).second.second.second);
			   if(m==io_map.end())std::cout << "Could not find Output signal: " << (*lr).second.second.second << " in hashtable." << std::endl;
			   else (*k).push_back(((*m).second.second));
		   }
	   }
	 }
	
	for(auto p = s.pi.cbegin();p!=s.pi.cend();p++)  // Iterator over all pattern 
	{
		// Build the expected output values for the primary outputs
		auto primary_output_data = (*p).multi_and_allclock_info.back().back().second.cbegin();
		for(auto i = s.signal_groups[1].signal.second.cbegin();i!= s.signal_groups[1].signal.second.cend();i++,primary_output_data++)
		{
			if((*i).find("test_so")==std::string::npos)
			{
				//lookup the name in the hashtable, if it does not exist create it
				bitset_map::iterator m =  po_expected_values.find((*i));
				if(m==po_expected_values.end())
				{
				//	std::cout << "Could not find primary output: " << (*i) << ", entry will be created. " << std::endl;
					m = po_expected_values.insert(bitset_map::value_type((*i),boost_bitset())).first;
				
				}
				// entry already exist in the hash table, just add the bit
				(*m).second.push_back((*primary_output_data)=='L'?false:true);
			}
		}

		// Set the data for the Pseudo Primary Outputs
		auto scan_chain_iterator = s.scan_chains.cbegin();
		for(auto testso = (*p).input_pattern_maping_test_so.cbegin();testso!=(*p).input_pattern_maping_test_so.cend();testso++,scan_chain_iterator++) // iterate over all test_si's
		{
		//	std::cout << s.scan_chains.size() << " " << (*p).input_pattern_maping_test_so.size() <<  std::endl;
			auto scan_chain_order_iterator = (*scan_chain_iterator).scan_structure.rbegin();
			for(auto test_so_data_it = (*testso).second.cbegin();test_so_data_it!=(*testso).second.cend();test_so_data_it++,scan_chain_order_iterator++)
			{
			//	std::cout << (*testso).second.size() << " " << (*scan_chain_iterator).scan_structure.size() <<  std::endl;
				auto hash_table_ff_lookup = ff_info_map.find((*scan_chain_order_iterator).ff_name);
				if(hash_table_ff_lookup==ff_info_map.end())std::cout << "Could not find: " << (*scan_chain_order_iterator).ff_name << " in Hashtable. " << std::endl;
				else 
				{
					//const unsigned int a = (*io_map.find((*hash_table_ff_lookup).second.first.second)).second.second;
					auto f = io_map.find((*hash_table_ff_lookup).second.first.second);
					if(f==io_map.end())
					{
						std::cout << "Colud not find entry: " << (*hash_table_ff_lookup).second.first.second << " in Hashtable." << std::endl;
					}
					const std::string  test = (*io_map.find((*hash_table_ff_lookup).second.first.second)).first;
					bitset_map::iterator m =  po_expected_values.find(test);
					if(m==po_expected_values.end())
					{
				//		std::cout << "Could not find PPO_0: " << test << ", entry will be created. " << std::endl;
						m = po_expected_values.insert(bitset_map::value_type(test,boost_bitset())).first;
					}
					(*m).second.push_back((*test_so_data_it)=='L'?false:true);
				}
			}
		}
		//std::cout << "PPO Loop finished. " << std::endl;


		// assign the values to the primary inputs
		auto pattern_it = (*p).multi_and_allclock_info[0][0].second.cbegin();
		auto pi_iter = s.signal_groups[0].signal.second.cbegin();
		//std::cout << (*p).multi_and_allclock_info[0][0].second.size() << " " << s.signal_groups[0].signal.second.size() << std::endl;
		for(/*auto pattern_it = (*p).multi_and_allclock_info[0][0].second.cbegin()*/;pattern_it!=(*p).multi_and_allclock_info[0][0].second.cend();pattern_it++,pi_iter++)
		{
			if((*pi_iter)!="clock" || (*pi_iter).find("test_si")==std::string::npos || (*pi_iter).find("test_se")==std::string::npos)
			{
				auto lookup_pi = io_map.find((*pi_iter));
				if(lookup_pi != io_map.end()) circuit_states1[(*lookup_pi).second.second].push_back((*pattern_it)=='0'?false:true);
			//	else std::cout << "Could not find primary input: " << (*pi_iter) << std::endl;
			}

		}

		auto scan_chain_iterator_lookup_it = scan_chain_struc_lookup.cbegin();
		for(auto testsi = (*p).input_pattern_maping_test_si.cbegin();testsi!=(*p).input_pattern_maping_test_si.cend();testsi++,scan_chain_iterator_lookup_it++) // iterate over all test_si's
		{
			auto scan_chain_iterator_hash_itertaor = (*scan_chain_iterator_lookup_it).cbegin();
			for(auto test_si_data_it = (*testsi).second.cbegin();test_si_data_it!=(*testsi).second.cend();test_si_data_it++,scan_chain_iterator_hash_itertaor++)
			{
				//std::cout << circuit_states1[(*scan_chain_iterator_hash_itertaor)] << std::endl;
				//std::cout << (*scan_chain_iterator_hash_itertaor) << " != " <<  circuit_states1.size() << std::endl;
		/*		if((*scan_chain_iterator_hash_itertaor) > circuit_states1.size())
				{
					std::cout << "Error" << std::endl;
					std::cout <<  (*scan_chain_iterator_lookup_it).size() << " " << (*testsi).second.size() << std::endl;
					std::cout << (*scan_chain_iterator_hash_itertaor) << " != " <<  circuit_states1.size() << std::endl;
				}*/
				circuit_states1[(*scan_chain_iterator_hash_itertaor)].push_back((*test_si_data_it)=='0'?false:true);
				//std::cout << circuit_states1[(*scan_chain_iterator_hash_itertaor)] << std::endl;
			}
		}
	}

	// At the end also the unload data must be set, it would be nice if this loop could be included in the big pattern loop
	auto testso_data_it = s.uload_data.test_so.cbegin();
	auto scan_chain_order_iterator = s.scan_chains.cbegin();
	for(;scan_chain_order_iterator!=s.scan_chains.cend();scan_chain_order_iterator++,testso_data_it++)
	{
		auto a = (*testso_data_it).second.cbegin();
		for(auto i = (*scan_chain_order_iterator).scan_structure.rbegin();i!=(*scan_chain_order_iterator).scan_structure.rend();i++,a++)
		{
			auto hash_table_ff_lookup = ff_info_map.find((*i).ff_name);
			if(hash_table_ff_lookup==ff_info_map.end())std::cout << "Could not find: " << (*i).ff_name << " in Hashtable. " << std::endl;
			else 
			{
				//unsigned int b = (*io_map.find((*hash_table_ff_lookup).second.first.second)).second.second;
				//circuit_states1[b].push_back((*a)=='L'?false:true);
				//const unsigned int ab = (*io_map.find((*hash_table_ff_lookup).second.first.second)).second.second;
				const std::string  test = (*io_map.find((*hash_table_ff_lookup).second.first.second)).first;
				bitset_map::iterator m =  po_expected_values.find(test);
				if(m==po_expected_values.end())
				{
			//		std::cout << "Could not find PPO_1: " << test << ", entry will be created. " << std::endl;
					m = po_expected_values.insert(bitset_map::value_type(test,boost_bitset())).first;
				}
				//std::cout << "Found Enrty. " << std::endl;
				(*m).second.push_back((*a)=='L'?false:true);
			}
		}
	}
}

// this function generated random patterns for the input nodes 
boost_bitset generate(boost_bitset::size_type n) // the parameter n determines how many patterns should be generated, this function returns a bitset, for this reason the OpenCL function needs to call the function to_ulong() on the returned bitset in order to convert it to an unsigned integer
{
	static std::random_device rd;                                     // set up the random number generator
	const auto bpb = boost_bitset::bits_per_block;                   
	static std::vector<boost_bitset::block_type> arr((n + bpb - 1) / bpb);
	static std::independent_bits_engine<std::minstd_rand, boost_bitset::bits_per_block, boost_bitset::block_width_type> engine(rd());
	std::generate(std::begin(arr), std::end(arr), std::ref(engine)); // generate random numbers for the whole bitset
	if (n < bpb) return boost_bitset(n, arr[0]); // if the n is smaller than a block of the bitset, just return the first block

	return boost_bitset(std::cbegin(arr), std::cend(arr)); // otherwise the whole bitset
}

// this function is used to generate random patterns for the OpenCL version of the algorithm
// is does not really generate pattern, it uses the "generate" function which was defined above
// it just puts the gernerated pattens in the right position, because each thread has its own random pattern
void generate_pattern_cl(const Graph &g,
	const unsigned int &final_number_of_nodes,
	const unsigned int &TOTAL_NUMBER_OF_ELEMENTS_CL,     // the TOTAL_NUMBER_OF_ELEMENTS_CL is the total number of nodes: NumberOfThreads * NumberOfNodesInGraph
	std::vector<cl_uint> &cl_pattern1,
	std::vector<cl_uint> &cl_pattern2)
{
	for (unsigned int topsorted_iterator = 0; topsorted_iterator < final_number_of_nodes;topsorted_iterator++)  // this loop is iterating over all nodes 
	{
		const auto in_degree = boost::in_degree(topsorted_iterator, g);  // get the in_degree of every node
		if (!in_degree)                                                  // if the node is an input generate two pattern for it
		{
			for (unsigned int i = 0; i<TOTAL_NUMBER_OF_ELEMENTS_CL; i += final_number_of_nodes)  // now do this for every thread, 
			{
				cl_pattern1[topsorted_iterator + i] = generate(32).to_ulong(); // 1st pattern call the function from above, 32 because every thread will simulate 32 pattern
				cl_pattern2[topsorted_iterator + i] = generate(32).to_ulong(); // 2nd pattern 
			}
			
		}
		else
		{
			cl_pattern1[topsorted_iterator] = 0;                               // if the node is not an input set it to zero, this is actually not necessary because it 
			cl_pattern2[topsorted_iterator] = 0;                               // will overwritten anymays during the calculation with the data 
		}
	}

}

// this function just copies the information about the in_degree and out_degree of every node to arrays
// this must be done because in OpenCL these funtions can not be accessed
void prepare_data_cl(const Graph &g,
	const unsigned int final_number_of_nodes,
	const std::vector<char> &gate_types_vec,
	std::vector<cl_uint> &cl_gate_top_ordered,
	std::vector<std::pair<cl_uint, cl_uint>> &cl_in_edge_top_ordered,
	std::vector<std::pair<cl_uint, cl_uint>> &cl_out_edges_top_ordered
	)
{
	for (unsigned int topsorted_iterator = 0; topsorted_iterator < final_number_of_nodes;topsorted_iterator++)
	{
		const auto in_degree = boost::in_degree(topsorted_iterator, g);
		if (in_degree)
		{

			cl_gate_top_ordered[topsorted_iterator] = ((cl_uint)gate_types_vec[topsorted_iterator]);
			
			for (auto eit = boost::in_edges(topsorted_iterator, g); eit.first != eit.second; ++eit.first)
			{
				if (in_degree == 1) // it's a not gate
					cl_in_edge_top_ordered[topsorted_iterator] = std::make_pair((*eit.first).m_source, (*eit.first).m_source);
				else			// all other gates
				{
					const auto a = (*eit.first).m_source;
					eit.first++;
					const auto b = (*eit.first).m_source;
					cl_in_edge_top_ordered[topsorted_iterator] = std::make_pair(a, b);
				}

			}
		}
		else
		{
			cl_gate_top_ordered[topsorted_iterator] = 0 ;                       // if node is an input gate_type is zero
			cl_in_edge_top_ordered[topsorted_iterator] = std::make_pair(0, 0);   // if node is an input save (-1,-1)
		}

		const auto out_degree = boost::out_degree(topsorted_iterator, g);
		if (out_degree)
		{
			auto out_edges = boost::out_edges(topsorted_iterator, g);
			if (out_degree == 1)
				cl_out_edges_top_ordered[topsorted_iterator] = std::make_pair((*out_edges.first).m_target, (*out_edges.first).m_target);
			else
			{
				const auto a = (*out_edges.first).m_target;
				out_edges.first++;
				const auto b = (*out_edges.first).m_target;
				cl_out_edges_top_ordered[topsorted_iterator] = std::make_pair(a, b);
			}
		}
		else // Node must be an final output
		{
			cl_out_edges_top_ordered[topsorted_iterator] = std::make_pair(0, 0);
		}
	}
}

// this function creates a complete toplogical order of all the data,
// that means the the gate type and the information about the edges 
// are also toplogical ordered 
void oder_data_cl(const Graph &g,
			   const std::vector<cl_uint> &cl_top_order,
			   const std::vector<char> &gate_types_vec,
			   std::vector<cl_uint> &cl_gate_top_ordered,
			   std::vector<std::pair<cl_uint, cl_uint>> &cl_in_edge_top_ordered,
			   std::vector<std::pair<cl_uint, cl_uint>> &cl_out_edges_top_ordered
			   )
{
	for (const auto topsorted_iterator : cl_top_order)
	{
		const auto in_degree = boost::in_degree(topsorted_iterator, g);
		if (in_degree)
		{

			cl_gate_top_ordered.emplace_back((cl_uint)gate_types_vec[topsorted_iterator]);
			for (auto eit = boost::in_edges(topsorted_iterator, g); eit.first != eit.second; ++eit.first)
			{
				if (in_degree == 1) // it's a not gate
					cl_in_edge_top_ordered.emplace_back(std::make_pair((*eit.first).m_source, (*eit.first).m_source));
				else			// all other gates
				{
					const auto a = (*eit.first).m_source;
					eit.first++;
					const auto b = (*eit.first).m_source;
					cl_in_edge_top_ordered.emplace_back(std::make_pair(a, b));
				}

			}
		}
		else
		{
			cl_gate_top_ordered.emplace_back(0);                       // if node is an input gate_type is zero
			cl_in_edge_top_ordered.emplace_back(std::make_pair(0, 0));   // if node is an input save (-1,-1)
		}

		const auto out_degree = boost::out_degree(topsorted_iterator, g);
		if (out_degree)
		{
			auto out_edges = boost::out_edges(topsorted_iterator, g);
			if (out_degree == 1)
				cl_out_edges_top_ordered.emplace_back(std::make_pair((*out_edges.first).m_target, (*out_edges.first).m_target));
			else
			{
				const auto a = (*out_edges.first).m_target;
				out_edges.first++;
				const auto b = (*out_edges.first).m_target;
				cl_out_edges_top_ordered.emplace_back(std::make_pair(a, b));
			}
		}
		else // Node must be an final output
		{
			cl_out_edges_top_ordered.emplace_back(std::make_pair(0, 0));
		}
	}
}


// this function simply "generates" 4-Bit random numbers, it was used during the debug process 
std::string generate4bit()
{
	static std::string numbers[] = { "0010","0001","0011","0101","0101","0110","0111","1001","1010","1000","1100","1101","1001","1110","1111","0010","0011","0011","0100","0101","0110","0111","1011","1010","1011","1100","1101","1000","1110","1111"};
	std::random_shuffle(std::begin(numbers),std::end(numbers));
	return numbers[4]; //7 4
}

// this functtion maps the gate type from a string to an integer
char gateTypeToInteger(const std::string &gate_type)
{
	if(gate_type.find("NAND")!=  std::string::npos)return 1;
	else if(gate_type.find("NOR")!=std::string::npos)return 2;
	else if(gate_type.find("XOR")!=std::string::npos)return 3;
	else if(gate_type.find("AND")!=std::string::npos)return 4;
	else if(gate_type.find("OR")!=std::string::npos)return 5;
	else if(gate_type.find("INV")!=std::string::npos)return 6;
	else return 0;
}

// this function is used to check for errors which may be caused by OpenCL functions
void checkOpenCLError(const cl_int  &error, const std::string &&buffername, const std::string &&operation)
{
	if (error != CL_SUCCESS)
	{
		std::cout << "ErrorCode: " << error << std::endl;
		std::cout << ">> " + buffername + " could not be " + operation << ": " << (error == CL_MEM_OBJECT_ALLOCATION_FAILURE ? " Allocation Failure":" ") <<   std::endl;
	}
}

// sample_graph_writer contains information which can be used the graph image which can be created by the program, the fontsize etc. can be set here
struct sample_graph_writer {
	sample_graph_writer(const unsigned int &threadid) :thread_id(threadid){};
	void operator()(std::ostream& out) const {
		out << "graph [label=\"Thread " << thread_id << "\", fontname=Arial, fontsize=30.0 ]" << std::endl;
		out << "node [ fontname=Arial, fontsize=25.0 ]" << std::endl;
	}
private:
	unsigned int thread_id;
};

// this function is used for logging information, this information is used to generate strings for each node
// these strings contain show the pattern and robustness etc.
void generate_log_info(const Graph &g,
	std::vector<unsigned int> top_sorted_graph,
	const std::vector<std::string> &vertex_names,
	const std::vector<std::string>  &gate_type,
	const std::vector<boost_bitset> &circuit_states1,
	const std::vector<boost_bitset> &circuit_states2,
	const std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &is_robust_or_non_robust_result,
	const std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &is_robust_result,
	const std::vector<boost_bitset> &is_robust_path_to_po, 
	const std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> &belongs_to_robust_path_final_result,
	std::vector<std::string> &node_info
	)
{
	std::for_each(std::begin(node_info), std::end(node_info), [](std::string &str){ str.reserve(150); });
	std::sort(std::begin(top_sorted_graph),std::end(top_sorted_graph));
	const std::size_t NUM_OF_BITES = 4;
		for (unsigned int i = 0;i<top_sorted_graph.size();i++) //(const auto i : top_sorted_graph)
		{
			const auto in_degree = boost::in_degree(i, g);

			node_info[i] = "[" + vertex_names[i] + ", " + gate_type[i] + "]";
			node_info[i] += "\nP1 : " + boost::to_string(circuit_states1[i]).substr(0, NUM_OF_BITES) + " P2: " + boost::to_string(circuit_states2[i]).substr(0, NUM_OF_BITES);

			if (in_degree)
			{
				node_info[i] += "\nisRobust: ";
				for (const auto j : is_robust_result[i])
					node_info[i] += vertex_names[j.first] + " " + boost::to_string(j.second).substr(0,NUM_OF_BITES) + " ";
			}

			if (in_degree)
			{
				node_info[i] += "\nIsRorNR: ";
				for (const auto j : is_robust_or_non_robust_result[i])
					node_info[i] += vertex_names[j.first] + " " + boost::to_string(j.second).substr(0, NUM_OF_BITES) + " ";
			}

			node_info[i] += "\nIsRPtoPo: " + boost::to_string(is_robust_path_to_po[i]).substr(0, NUM_OF_BITES);

			node_info[i] += "\nBToRPath: ";
			for (const auto j : belongs_to_robust_path_final_result[i])
				node_info[i] += vertex_names[j.first] + " " + boost::to_string(j.second).substr(0, NUM_OF_BITES) + " ";
		}
}



#endif