/*
	This is the main file of the program. This file contains the flow of the program.
	It calls mostly functions (file parsing, path classification etc.) from the other files.
	This fle may look big but it just contains the flow and setting up for the functions. 
	Especially the code for setting up OpenCL is huge, a detailed explanaition about how to set up OpenCL can
	be found in the OpenCL book in the lab. 
*/

#include <boost/config/warning_disable.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/algorithm/string.hpp>
#define BOOST_DATE_TIME_NO_LIB
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/tuple/tuple.hpp>
#include <CL/cl.hpp>
#include <unordered_map>
#include <string>
#include <utility>
#include <fstream>
#include <cassert>
#include <iostream>  
#include <vector>
#include <future>
#include <algorithm>
#include <iterator>
#include <functional>
#include <exception>

#include "file_parsers.h"
#include "graph_structure.h"
#include "graph_creator.h"
#include "hash_maps.h"
#include "helper_functions.h"
#include "pathclassification.h"
#include "logicsimulation_cpu.h"

const bool CPU_LOGIC_AND_PATH_SIM = false;
const bool OPENCL_SIMULATION = true;
const bool OPENCL_GPU_PATH_CLASSIFICATION = false; // if this is set to false and OPENCL_SIMULATION is true then OpenCL CPU is used 
const bool STILEFILE_CHECK = false;
const bool LOGGING = false;
const bool OPENCL_LOGGING = false;

std::size_t NUMBER_OF_THREADS = 512;
const std::size_t WORK_GROUP_SIZE = 1;

boost_bitset::size_type bppk = 32768; // 16384

int main(int argc, char* argv[])
{
	boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::universal_time();
	if (argc < 3)
	{
		std::cout << "Usage: <Netlist>.vg <Stilfile>.stil number_of_threads" << std::endl;
		return 0;

	}
	const std::string filename(argv[1]);  //leon3mp_scan
	const std::string stilfilename(argv[2]); // leon3mp_tdf.stil
	NUMBER_OF_THREADS = std::stoi(argv[3]);
	
	if (CPU_LOGIC_AND_PATH_SIM)
	{
		bppk = NUMBER_OF_THREADS * 32;
	}
	std::cout << "Numer of pattern: " << bppk << std::endl;

	if (!std::ifstream(filename.c_str()).good() || !std::ifstream(stilfilename.c_str()).good())
	{
		std::cout << "One or both of the files do not exist. " << std::endl;
		return 0;
	}

	std::future<stilfile_datastruct> s;
	if (STILEFILE_CHECK)
		s = std::async([&stilfilename]() { return parse_stilfile(stilfilename); });

	std::future<circuit> c = std::async(std::launch::async,[&filename]() { return parse_netlist(filename); });


	my_map io_map;
	dff_map ff_info_map;
	unsigned int counter = 0;

	std::vector<char>gate_types_vec;
	std::vector<std::string>vertex_names;
	std::vector<std::string>gate_type;
	std::vector<std::pair<bool, unsigned int>> constant_input_or_wires;
	// Create the graph edges
	unsigned int final_number_of_nodes;
	Graph g = create_graph(c.get(), vertex_names, gate_type, gate_types_vec, constant_input_or_wires, final_number_of_nodes,io_map,ff_info_map);

	std::vector<unsigned int> top_sorted_graph;
	top_sorted_graph.reserve(final_number_of_nodes);
	boost::topological_sort(g, std::back_inserter(top_sorted_graph)); 
	std::reverse(std::begin(top_sorted_graph), std::end(top_sorted_graph));

	Graph transponse_graph;
	boost::transpose_graph(g, transponse_graph);


	if (CPU_LOGIC_AND_PATH_SIM)
	{
		std::vector<boost_bitset>circuit_states1(final_number_of_nodes);
		std::vector<boost_bitset>circuit_states2(final_number_of_nodes);
		std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> is_robust_or_non_robust_result;
		std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> is_robust_result;
		std::vector<std::vector<std::pair<unsigned int, boost_bitset>>> belongs_to_robust_path_final_result;
		std::vector<boost_bitset> is_robust_path_to_po;

		stilfile_datastruct s1;
		if (STILEFILE_CHECK)
		{
			s1 = s.get();
			bppk = s1.pi.size();
		}



		//boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::universal_time();
		if (STILEFILE_CHECK)
		{
			bitset_map po_expected_values;
			stilfile_check(s1, io_map, ff_info_map, circuit_states1, po_expected_values);
			circuit_states2 = circuit_states1;
			logic_simlation_cpu(g, top_sorted_graph, gate_types_vec, bppk, circuit_states1, circuit_states2);
			std::cout << "Logic Simulation 1 finished " << std::endl;
			for (const auto i : ff_info_map)
			{
				auto hash_table_lookup_ppo = i.second.first.second;
				auto hash_table_lookup_ppi = i.second.second.second;
				auto lookup_result1 = io_map.find(hash_table_lookup_ppo);
				auto lookup_result2 = io_map.find(hash_table_lookup_ppi);
				if (lookup_result1 == io_map.end() || lookup_result2 == io_map.end()) std::cout << "Could not find: " << hash_table_lookup_ppo << std::endl;
				else
				{
					circuit_states1[(*lookup_result2).second.second] = circuit_states1[(*lookup_result1).second.second];
					circuit_states2[(*lookup_result2).second.second] = circuit_states2[(*lookup_result1).second.second];
				}
			}

			logic_simlation_cpu(g, top_sorted_graph, gate_types_vec, bppk, circuit_states1, circuit_states2);

			for (const auto i : po_expected_values)
			{
				const auto res_lookup = io_map.find(i.first);
				if (res_lookup == io_map.end())
					std::cout << "Could not find entry in table. " << std::endl;
				else
				{
					if (circuit_states1[(*res_lookup).second.second] != (i).second)
					{
						std::cout << "MISMATCH for " << i.first << std::endl; // " expected value: " << (*i).second << " simulated value: " <<  circuit_states1[(*res_lookup).second.second] << std::endl;
					}
				}

			}
		}
		else
		{

			for (unsigned int i = 0; i < final_number_of_nodes;i++) //const auto i : top_sorted_graph
			{
				if (!boost::in_degree(i, g))
				{
					circuit_states1[i] = generate(bppk);
					circuit_states2[i] = generate(bppk);
				}
			}
			for (const auto i : constant_input_or_wires)
			{
				if (i.first)
				{
					std::cout << "Setting Wire to 1" << std::endl;
					circuit_states1[i.second] = circuit_states2[i.second] = boost_bitset(bppk).set();
				}
				else
				{
					std::cout << "Setting Wire to 0" << std::endl;
					circuit_states1[i.second] = circuit_states2[i.second] = boost_bitset(bppk);
				}
			}
			boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::universal_time();
			std::cout << "Preprocessing finished after:" << (t2 - t1).total_milliseconds() << " Milliseconds preprocessing" << std::endl;
			//boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::universal_time();
			//logic_simlation_cpu(g, top_sorted_graph, gate_types_vec, bppk, circuit_states1, circuit_states2);	// The Logic Simulation
			//logic_simlation_cpu_top(g, final_number_of_nodes, gate_types_vec, bppk, circuit_states1, circuit_states2);
			
			path_classification_with_logic(g, transponse_graph, circuit_states1, circuit_states2, gate_types_vec, final_number_of_nodes,
				bppk, vertex_names, is_robust_or_non_robust_result, is_robust_result, is_robust_path_to_po, belongs_to_robust_path_final_result);

			//boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::universal_time();
			std::cout << "Simulation finished after:" << (t2 - t1).total_milliseconds() << " Milliseconds preprocessing" << std::endl;

			if (LOGGING)
			{
				std::vector<std::string> node_info(final_number_of_nodes);
				generate_log_info(g,
					top_sorted_graph,
					vertex_names,
					gate_type,
					circuit_states1,
					circuit_states2,
					is_robust_or_non_robust_result,
					is_robust_result,
					is_robust_path_to_po,
					belongs_to_robust_path_final_result,
					node_info);
				std::ofstream output("Graph_CPU.dot");
				boost::write_graphviz(output, g, boost::make_label_writer(node_info.data()));
			}
		}
	}
	boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::universal_time();
	//std::cout << "Simulation finished:" << (t2 - t1).total_microseconds() << " Microseconds" << std::endl;

	// From here the OpenCL program starts **************************************************************************************************************
	if (OPENCL_SIMULATION)
	{
		std::cout << "Number of threads:  " << NUMBER_OF_THREADS << std::endl;
		const auto TOTAL_NUMBER_OF_ELEMENTS_CL = NUMBER_OF_THREADS * final_number_of_nodes;

	//	std::vector<cl_uint>cl_top_order(std::cbegin(top_sorted_graph), std::cend(top_sorted_graph));
	//	std::vector<cl_uint>cl_revtop_order(std::cbegin(top_sorted_graph), std::cend(top_sorted_graph));
		std::vector<cl_uint>cl_gates(final_number_of_nodes);


		std::vector<cl_uint>cl_pattern1(TOTAL_NUMBER_OF_ELEMENTS_CL);   // not toplogical sorted
		std::vector<cl_uint>cl_pattern2(TOTAL_NUMBER_OF_ELEMENTS_CL);   // not toplogical sorted

		std::vector<std::pair<cl_uint, cl_uint>> cl_in_edges(final_number_of_nodes);
		std::vector<std::pair<cl_uint, cl_uint>> cl_out_edges(final_number_of_nodes);

		// constant wires must be set before the simulation starts for every thread
	


		prepare_data_cl(g, final_number_of_nodes, gate_types_vec, cl_gates, cl_in_edges, cl_out_edges);
		generate_pattern_cl(g, final_number_of_nodes, TOTAL_NUMBER_OF_ELEMENTS_CL, cl_pattern1, cl_pattern2);

		for (const auto i : constant_input_or_wires)
			for (std::size_t k = 0; k < TOTAL_NUMBER_OF_ELEMENTS_CL; k += final_number_of_nodes)
			{
				//std::cout << "Setting constant values" << std::endl;
				cl_pattern1[i.second + k] = cl_pattern2[i.second + k] = i.first ? ~0 : 0;
			}

		std::cout << "Memory allocation starts." << TOTAL_NUMBER_OF_ELEMENTS_CL << std::endl;
		boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::universal_time();
		std::vector<cl_uint> cl_is_robust1(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_is_robust2(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_is_robust_or_non_robust1(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_is_robust_or_non_robust2(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_is_robust_path_to_po(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_belongs_to_robust_path_final_result1(TOTAL_NUMBER_OF_ELEMENTS_CL);
		std::vector<cl_uint> cl_belongs_to_robust_path_final_result2(TOTAL_NUMBER_OF_ELEMENTS_CL);
		boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::universal_time();
		std::cout << "Memory Allocation finished after:" << (t2 - t1).total_milliseconds() << " Milliseconds preprocessing" << std::endl;

		/*
		
			From here the OpenCL code starts
		*/
		std::vector<cl::Platform>platforms;
		cl::Platform::get(&platforms);
		for (auto i : platforms)
			std::cout << i.getInfo<CL_PLATFORM_VENDOR>() << std::endl;

		cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[1])(), 0 };
		cl::Context context(CL_DEVICE_TYPE_ALL, cps);
		std::vector<cl::Device>devices = context.getInfo<CL_CONTEXT_DEVICES>();
		std::cout << "Enumerating all OpenCL devices: " << std::endl;
		for (const auto dev_it : devices)
			std::cout << "   " << dev_it.getInfo<CL_DEVICE_NAME>()  <<  std::endl;

		cl_int err;
		cl::CommandQueue queue;
		if (devices.size() >= 2)
		{
			std::cout << "Number of Devcices > 2 " << std::endl;
			if (OPENCL_GPU_PATH_CLASSIFICATION)
				queue = cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE,&err);
			else
				queue = cl::CommandQueue(context, devices[1], CL_QUEUE_PROFILING_ENABLE, &err);
		}
		else queue = cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &err);
		checkOpenCLError(err, "CommandQueue", "created");


		cl_mem_flags cbuffer_flags = CL_MEM_READ_ONLY;
		cl_mem_flags buffer_flags  = CL_MEM_READ_WRITE;
		if (!OPENCL_GPU_PATH_CLASSIFICATION)
		{
			cbuffer_flags |= CL_MEM_USE_HOST_PTR;
			buffer_flags |= CL_MEM_USE_HOST_PTR;
		}
		else
		{
			cbuffer_flags |= CL_MEM_COPY_HOST_PTR;
			buffer_flags |= CL_MEM_COPY_HOST_PTR;
		}

		// Create the buffers 
	//	cl::Buffer top_order_buffer = cl::Buffer(context, cbuffer_flags, cl_top_order.size()*sizeof(*cl_top_order.data()), cl_top_order.data(), &err);           checkOpenCLError(err, "top_order_buffer", "created");
	//	cl::Buffer revtop_order_buffer = cl::Buffer(context, cbuffer_flags, cl_revtop_order.size()*sizeof(*cl_revtop_order.data()), cl_revtop_order.data(), &err);   checkOpenCLError(err, "top_reverse_order_buffer", "created");

		cl::Buffer gate_type_buffer = cl::Buffer(context, cbuffer_flags, cl_gates.size()*sizeof(*cl_gates.data()), cl_gates.data(), &err);        checkOpenCLError(err, "gate_type_buffer", "created");
		cl::Buffer in_edges_buffer = cl::Buffer(context, cbuffer_flags, cl_in_edges.size()*sizeof(*cl_in_edges.data()), cl_in_edges.data(), &err);       checkOpenCLError(err, "inedges_buffer", "created");
		cl::Buffer out_edges_buffer = cl::Buffer(context, cbuffer_flags, cl_out_edges.size()*sizeof(*cl_out_edges.data()), cl_out_edges.data(), &err);       checkOpenCLError(err, "outedges_buffer", "created");

		cl::Buffer pattern_buffer1 = cl::Buffer(context, buffer_flags, cl_pattern1.size() * sizeof(*cl_pattern1.data()), cl_pattern1.data(), &err); checkOpenCLError(err, "pattern_buffer1", "created");
		cl::Buffer pattern_buffer2 = cl::Buffer(context, buffer_flags, cl_pattern2.size() * sizeof(*cl_pattern2.data()), cl_pattern2.data(), &err); checkOpenCLError(err, "pattern_buffer2", "created");
		cl::Buffer robust_buffer1 = cl::Buffer(context, buffer_flags, cl_is_robust1.size() * sizeof(*cl_is_robust1.data()), cl_is_robust1.data(), &err); checkOpenCLError(err, "robust_buffer1", "created");
		cl::Buffer robust_buffer2 = cl::Buffer(context, buffer_flags, cl_is_robust2.size() * sizeof(*cl_is_robust2.data()), cl_is_robust2.data(), &err); checkOpenCLError(err, "robust_buffer2", "created");
		cl::Buffer nonrobust_buffer1 = cl::Buffer(context, buffer_flags, cl_is_robust_or_non_robust1.size() * sizeof(*cl_is_robust_or_non_robust1.data()), cl_is_robust_or_non_robust1.data(), &err); checkOpenCLError(err, "is_robust_or_non_robuts_buffer1", "created");
		cl::Buffer nonrobust_buffer2 = cl::Buffer(context, buffer_flags, cl_is_robust_or_non_robust2.size() * sizeof(*cl_is_robust_or_non_robust2.data()), cl_is_robust_or_non_robust2.data(), &err); checkOpenCLError(err, "is_robust_or_non_robuts_buffer2", "created");
		cl::Buffer robustpathtopo_buffer = cl::Buffer(context, buffer_flags, cl_is_robust_path_to_po.size() * sizeof(*cl_is_robust_path_to_po.data()), cl_is_robust_path_to_po.data(), &err); checkOpenCLError(err, "is_robust_path_to_po", "created");
		cl::Buffer belongsrpath_buffer1 = cl::Buffer(context, buffer_flags, cl_belongs_to_robust_path_final_result1.size() * sizeof(*cl_belongs_to_robust_path_final_result1.data()), cl_belongs_to_robust_path_final_result1.data(), &err); checkOpenCLError(err, "bto_robust_path_buffer", "created");
		cl::Buffer belongsrpath_buffer2 = cl::Buffer(context, buffer_flags, cl_belongs_to_robust_path_final_result2.size() * sizeof(*cl_belongs_to_robust_path_final_result2.data()), cl_belongs_to_robust_path_final_result2.data(), &err); checkOpenCLError(err, "bto_robust_path_buffer", "created");

		if (OPENCL_GPU_PATH_CLASSIFICATION)
		{
	//		err = queue.enqueueWriteBuffer(top_order_buffer, CL_FALSE, 0, cl_top_order.size()*sizeof(*cl_top_order.data()), &cl_top_order.front()); checkOpenCLError(err, "top_order_buffer", "transfer");
	//		err = queue.enqueueWriteBuffer(revtop_order_buffer, CL_FALSE, 0, cl_revtop_order.size()*sizeof(*cl_revtop_order.data()), &cl_revtop_order.front()); checkOpenCLError(err, "rev_top_order_buffer", "transfer");
			err = queue.enqueueWriteBuffer(gate_type_buffer, CL_FALSE, 0, cl_gates.size()*sizeof(*cl_gates.data()), &cl_gates.front());  checkOpenCLError(err, "gate_type_buffer", "transfer");
			err = queue.enqueueWriteBuffer(in_edges_buffer, CL_FALSE, 0, cl_in_edges.size()*sizeof(*cl_in_edges.data()), &cl_in_edges.front()); checkOpenCLError(err, "in_edges_top_buffer", "transfer");
			err = queue.enqueueWriteBuffer(out_edges_buffer, CL_FALSE, 0, cl_out_edges.size()*sizeof(*cl_out_edges.data()), &cl_out_edges.front()); checkOpenCLError(err, "out_edge_rev__buffer", "transfer");



			err = queue.enqueueWriteBuffer(pattern_buffer1, CL_FALSE, 0, cl_pattern1.size() * sizeof(*cl_pattern1.data()), &cl_pattern1.front()); checkOpenCLError(err, "pattern_buffer1", "transfer");
			err = queue.enqueueWriteBuffer(pattern_buffer2, CL_FALSE, 0, cl_pattern2.size() * sizeof(*cl_pattern2.data()), &cl_pattern2.front()); checkOpenCLError(err, "pattern_buffer2", "transfer");

			err = queue.enqueueWriteBuffer(robust_buffer1, CL_FALSE, 0, cl_is_robust1.size() * sizeof(*cl_is_robust1.data()), &cl_is_robust1.front()); checkOpenCLError(err, "robust_buffer1", "transfer");
			err = queue.enqueueWriteBuffer(robust_buffer2, CL_FALSE, 0, cl_is_robust2.size() * sizeof(*cl_is_robust2.data()), &cl_is_robust2.front()); checkOpenCLError(err, "robust_buffer2", "transfer");

			err = queue.enqueueWriteBuffer(nonrobust_buffer1, CL_FALSE, 0, cl_is_robust_or_non_robust1.size() *sizeof(*cl_is_robust_or_non_robust1.data()), &cl_is_robust_or_non_robust1.front()); checkOpenCLError(err, "cl_is_robust_or_non_robust_buffer1", "transfer");
			err = queue.enqueueWriteBuffer(nonrobust_buffer2, CL_FALSE, 0, cl_is_robust_or_non_robust2.size() *sizeof(*cl_is_robust_or_non_robust2.data()), &cl_is_robust_or_non_robust2.front()); checkOpenCLError(err, "cl_is_robust_or_non_robust_buffer2", "transfer");

			err = queue.enqueueWriteBuffer(robustpathtopo_buffer, CL_FALSE, 0, cl_is_robust_path_to_po.size() *sizeof(*cl_is_robust_path_to_po.data()), &cl_is_robust_path_to_po.front()); checkOpenCLError(err, "robust_path_to_po_buffer", "transfer");

			err = queue.enqueueWriteBuffer(belongsrpath_buffer1, CL_FALSE, 0, cl_belongs_to_robust_path_final_result1.size() *sizeof(*cl_belongs_to_robust_path_final_result1.data()), &cl_belongs_to_robust_path_final_result1.front());
			checkOpenCLError(err, "belongsrpath_buffer1", "transfer");
			err = queue.enqueueWriteBuffer(belongsrpath_buffer2, CL_FALSE, 0, cl_belongs_to_robust_path_final_result2.size() *sizeof(*cl_belongs_to_robust_path_final_result2.data()), &cl_belongs_to_robust_path_final_result2.front());
			checkOpenCLError(err, "belongsrpath_buffer2", "transfer");
		}


		std::ifstream sourceFile("PathDelayFault.cl");
		std::string sourceCode(std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));

		cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length() + 1));
		cl::Program program = cl::Program(context, source);
		// Build program for these specific devices
		program.build(devices, std::string("-DSIZE=" + std::to_string(final_number_of_nodes)).c_str());

		// Make kernel
		const std::string kernel_name = "logicsimulationv2";
		std::cout << "Executing Kernel: " << kernel_name << std::endl;
		cl::Kernel kernel(program, kernel_name.c_str(), &err);
		for (const auto dev : devices)
		{
			std::cout << "I like work-group sizes of: " << kernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(dev) << std::endl;
		}
		
		checkOpenCLError(err, "Kernel", "be compiled");

		const std::vector<cl::Buffer*>set_args = { /*&top_order_buffer,
			&revtop_order_buffer,*/
			&gate_type_buffer,
			&in_edges_buffer,
			&out_edges_buffer,
			&pattern_buffer1,
			&pattern_buffer2,
			&robust_buffer1,
			&robust_buffer2,
			&nonrobust_buffer1,
			&nonrobust_buffer2,
			&robustpathtopo_buffer,
			&belongsrpath_buffer1,
			&belongsrpath_buffer2 };

	//	kernel.setArg(0, top_sorted_graph.size());
		for (cl_uint i = 0; i < set_args.size(); i++)
			kernel.setArg(i, *set_args[i]);


		cl::NDRange global(NUMBER_OF_THREADS);
		cl::NDRange local(WORK_GROUP_SIZE);
		//cl::NullRange 
		cl::Event event;
		err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, global,local , nullptr, &event);
		checkOpenCLError(err, "kernel", "executed");
			
		cl_ulong  queued = 0;
		cl_ulong  submit = 0;
		cl_ulong  start = 0;
		cl_ulong  end = 0;
		event.wait();
		event.getProfilingInfo<cl_ulong>(CL_PROFILING_COMMAND_QUEUED,&queued);
		event.getProfilingInfo<cl_ulong>(CL_PROFILING_COMMAND_SUBMIT, &submit);
		event.getProfilingInfo<cl_ulong>(CL_PROFILING_COMMAND_START, &start);
		event.getProfilingInfo<cl_ulong>(CL_PROFILING_COMMAND_END, &end);
		std::ofstream profile_data("b19_cpu.csv", std::ios::app);
		//std::cout << queued << " " << submit << " " << submit-queued << start << " " << start-submit << " " << end << " " << end-start << std::endl;
		profile_data << (end - start) << std::endl;
		

		if (OPENCL_GPU_PATH_CLASSIFICATION)
		{
			queue.enqueueReadBuffer(pattern_buffer1, CL_TRUE, 0, cl_pattern1.size() * sizeof(*cl_pattern1.data()), &cl_pattern1.front());
			queue.enqueueReadBuffer(pattern_buffer2, CL_TRUE, 0, cl_pattern2.size() * sizeof(*cl_pattern2.data()), &cl_pattern2.front());
			queue.enqueueReadBuffer(robust_buffer1, CL_TRUE, 0, cl_is_robust1.size() * sizeof(*cl_is_robust1.data()), &cl_is_robust1.front());
			queue.enqueueReadBuffer(robust_buffer2, CL_TRUE, 0, cl_is_robust2.size() * sizeof(*cl_is_robust2.data()), &cl_is_robust2.front());
			queue.enqueueReadBuffer(nonrobust_buffer1, CL_TRUE, 0, cl_is_robust_or_non_robust1.size() * sizeof(*cl_is_robust_or_non_robust1.data()), &cl_is_robust_or_non_robust1.front());
			queue.enqueueReadBuffer(nonrobust_buffer2, CL_TRUE, 0, cl_is_robust_or_non_robust2.size() * sizeof(*cl_is_robust_or_non_robust2.data()), &cl_is_robust_or_non_robust2.front());
			queue.enqueueReadBuffer(robustpathtopo_buffer, CL_TRUE, 0, cl_is_robust_path_to_po.size() * sizeof(*cl_is_robust_path_to_po.data()), &cl_is_robust_path_to_po.front());
			queue.enqueueReadBuffer(belongsrpath_buffer1, CL_TRUE, 0, cl_belongs_to_robust_path_final_result1.size() * sizeof(*cl_belongs_to_robust_path_final_result1.data()), &cl_belongs_to_robust_path_final_result1.front());
			queue.enqueueReadBuffer(belongsrpath_buffer2, CL_TRUE, 0, cl_belongs_to_robust_path_final_result2.size() * sizeof(*cl_belongs_to_robust_path_final_result2.data()), &cl_belongs_to_robust_path_final_result2.front()); 
		}

		if (OPENCL_LOGGING)
		{
			const std::size_t OUTPUT_BITS = 4;
			std::vector<std::string>gpu_gate_type_and_state(final_number_of_nodes);
			std::for_each(std::begin(gpu_gate_type_and_state), std::end(gpu_gate_type_and_state), [](std::string &str){str.reserve(170); });

			const std::string dot_path = "\"C:\\Users\\user\\Desktop\\release\\bin\\dot.exe\"";
			const std::string filename_wihtout_dot = std::string(filename.begin(), filename.begin() + filename.find_first_of('.'));
			const std::string graph_dic = filename_wihtout_dot + "_Images";

			std::vector<std::string> dot_files = { "@echo off", "set dotpath=" + dot_path, "mkdir " + graph_dic };

			std::cout << "Writing GPU-Graphs to file: ";
			for (auto k = 0; k < NUMBER_OF_THREADS; k++)
			{
				const unsigned int print_thread_data_nr = k * final_number_of_nodes;
				//for (const auto i : cl_top_order)
				for (unsigned int i = 0; i < final_number_of_nodes;i++)
				{
					bool input = false;
					bool output = false;
					gpu_gate_type_and_state[i] = gate_type[i] + " " + vertex_names[i];

					if (cl_in_edges[i].first != cl_in_edges[i].second) gpu_gate_type_and_state[i] += "\nInEdges: " + vertex_names[cl_in_edges[i].first] + "  " + vertex_names[cl_in_edges[i].second];
					else if ((cl_in_edges[i].first == cl_in_edges[i].second) && (cl_in_edges[i].first != 0 && cl_in_edges[i].second != 0))
					{

						gpu_gate_type_and_state[i] += "\nInEdges: " + vertex_names[cl_in_edges[i].first];
					}
					else input = true;



					if (cl_out_edges[i].first != cl_out_edges[i].second)
						gpu_gate_type_and_state[i] += "\nOutEdges: " + vertex_names[cl_out_edges[i].first] + "  " + vertex_names[cl_out_edges[i].second];
					else if ((cl_out_edges[i].first == cl_out_edges[i].second) && cl_out_edges[i].first != 0 && cl_out_edges[i].second != 0)
					{
						gpu_gate_type_and_state[i] += "\nOutEdges: " + vertex_names[cl_out_edges[i].first];
					}
					else output = true;

					gpu_gate_type_and_state[i] += "\nP1: " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_pattern1[i + print_thread_data_nr]));
					gpu_gate_type_and_state[i] += " P2: " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_pattern2[i + print_thread_data_nr]));

					if (!input)
					{
						gpu_gate_type_and_state[i] += "\nisR: " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_is_robust1[i + print_thread_data_nr])) + " " +
							boost::to_string(boost_bitset(OUTPUT_BITS, cl_is_robust2[i + print_thread_data_nr]));
						//"\nisRorNR: " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_is_robust_or_non_robust1[c_node + print_thread_data_nr])) + "  " +
						//boost::to_string(boost_bitset(OUTPUT_BITS, cl_is_robust_or_non_robust2[c_node + print_thread_data_nr])) +
					}


					gpu_gate_type_and_state[i] += "\nisRobustPathTo: " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_is_robust_path_to_po[i + print_thread_data_nr]));

					gpu_gate_type_and_state[i] += "\nBToRPath: ";/* +boost::to_string(boost_bitset(OUTPUT_BITS, cl_pattern1[i + print_thread_data_nr]));
					gpu_gate_type_and_state[i] += " " + boost::to_string(boost_bitset(OUTPUT_BITS, cl_pattern2[i + print_thread_data_nr]));*/
					if (cl_out_edges[i].first != cl_out_edges[i].second)
					{
						gpu_gate_type_and_state[i] += boost::to_string(boost_bitset(OUTPUT_BITS, cl_belongs_to_robust_path_final_result1[i + print_thread_data_nr])) + " " +
						boost::to_string(boost_bitset(OUTPUT_BITS, cl_belongs_to_robust_path_final_result2[i + print_thread_data_nr]));
					}
							else
							gpu_gate_type_and_state[i] += boost::to_string(boost_bitset(OUTPUT_BITS, cl_belongs_to_robust_path_final_result1[i + print_thread_data_nr]));
				}
				std::string outputfile = "Graph_thread" + std::to_string(k);
				std::ofstream output(std::string(outputfile + ".dot").c_str());
				dot_files.emplace_back("%dotpath%  -Tpng " + outputfile + ".dot > " + outputfile + ".png");
				dot_files.emplace_back("move " + outputfile + ".png " + graph_dic);
				dot_files.emplace_back("del " + outputfile + ".dot ");
				//boost::write_graphviz(output, g, );
				std::cout << k << " ";
				boost::write_graphviz(output, g, boost::make_label_writer(&gpu_gate_type_and_state[0]), boost::default_writer(), sample_graph_writer(k));
				output.close();

			}
			std::copy(std::cbegin(dot_files), cend(dot_files), std::ostream_iterator<std::string>(std::ofstream(std::string(filename_wihtout_dot + "_image_script.bat").c_str()), "\n"));
			std::cout << std::endl;
		}
	}
	std::cout << "Programmende. " << std::endl;
	//std::cin.get();
	return 0;
}


