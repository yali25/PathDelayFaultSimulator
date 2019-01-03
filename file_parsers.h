/*
  This file is used for parsing the data from the Netlist-file into a datastructure which
  makes it comfortable to work with with. 
  This file contains two parsers, Netlist-Parser (vgdata)  and Stilefile-Parser (stilfileparser).
  This part of the program makes heavy usage of the Boost Spirit Library. Please check the 
  boost Spirit documentation for details.
*/

#ifndef FILE_PARSER
#define FILE_PARSER

#include <string>
#include <vector>
#include <boost/spirit/include/qi.hpp>  
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/dynamic_bitset.hpp>

/*
From here the Netlist-Parser part starts.
*/
// This typedef is used to represent the inputs and outputs of a gate.
// For Example the gate: NAND2X1 U16 ( .IN1(G6), .IN2(n3), .QN(n6) )
// has two inputs and one ouput, so this datastructure will hold the data:
//  { IN1(G6), IN2(n3), QN(n6) }
// the std::pair datastructure is a built in type of C++, its data can be accessed
// by "first" and "second"
// In the example above three pairs are inserted into the vector { (IN1,G6),(IN2,n3),(QN,n6) } 
typedef std::vector<std::pair<std::string, std::string>> vec_type;

// This data structure is used to represent a gate.
// a gate has a type e.g. AND, OR etc.
// a gate has also name e.g. U16
// and each gate has list of inputs ans outputs of the type "vec_type" which was defined above "ipc" stands for "input pin configuration"
struct gate
{
	std::string gate_type;
	std::string gate_name;
	vec_type ipc;
};
// this macro (BOOST_FUSION_ADAPT_STRUCT) must be called to tell the boost Spirit library how the data-strcture looks like/
// first parameter is the name of the struct, and then all datafields and there types must be written
BOOST_FUSION_ADAPT_STRUCT(gate, (std::string, gate_type)
	(std::string, gate_name)
	(vec_type, ipc)
	)

	// the circuit_io_object data structure is used to represent "input" and "output" elements of
	// a netlist file. "type" saves if it is an "input" or "ouput" and then a list of names is saved
	// in the vector io_elements, this information is actually never used in the program but maybe
	// it is useful for future research projects
struct circuit_io_object
{
	std::string type;
	std::vector<std::string> io_elements;
};
//  also call the macro to tell Boost how the data structure looks
BOOST_FUSION_ADAPT_STRUCT(circuit_io_object, (std::string, type)
	(std::vector<std::string>, io_elements)
	)

	// this type was defined for saving the wires of the netlist file, but actually it never used,
	// maybe it may also be useful for the future
	typedef std::vector<std::string> wire_vec;

// The circuit data-structure represents the final data of a whole netlist file.
// All previously defined data structures are used here. A circuit has a name e.g. "s27"
// and a list of inputs e.g. s27 (g1, j5, h5....)
// and a list of Input-Output Elements 
// and the most important information all the gates of the netlist
// NOTE: there is no list of wires. If the list of wires is needed 
// a new element new list of type  wire_vec must be added
struct circuit
{
	std::string                    circuit_name;
	std::vector<std::string>       circuit_inputs;
	std::vector<circuit_io_object> circuit_io_elements;
	std::vector<gate>              circuit_gates;

};
// use again the macro to tell boost how the data structure looks
BOOST_FUSION_ADAPT_STRUCT(circuit, (std::string, circuit_name)
	(std::vector<std::string>, circuit_inputs)
	(std::vector<circuit_io_object>, circuit_io_elements)
	(std::vector<gate>, circuit_gates)
	)
	// From here the actuall parser starts, for understanding how this works it is useful
	// to understand the Backus-Naur_Form  (http://en.wikipedia.org/wiki/Backus-Naur_Form)
	// I will explain hoe the rule identifier %=  qi::lexeme[ (ascii::alnum >> -qi::char_('\'') >> *qi::char_("a-zA-Z_0-9") >> -(qi::char_('[') >> *ascii::digit >> qi::char_(']'))) ] works
	// qi::lexeme says that as this level no spaces should be skipped, becaued if they would be skipped a wrong netlist could be parsed.
	// ascii::alnum says that the identifier can start with a number or a character, then a ' may follow, this is because e.g. 1'G1 is valid in the netlist file, it means G1 is constant 1
	// *qi::char_("a-zA-Z_0-9") says that now an abritraly long string with mixed characters and numbers can follow, also *ascii::alnum could be used here but the name could also contain _ and that character is not covered by ascii::alnum
	// -(qi::char_('[') >> *ascii::digit >> qi::char_(']')) says that there may an open bracket followed with an integer and a closing bracket because b1[12] is a valid name
	template <typename Iterator, typename Skipper>
struct vgdata : boost::spirit::qi::grammar<Iterator, circuit(), Skipper>
{
	vgdata() : vgdata::base_type(object)
	{
		using namespace boost::spirit;
		identifier %= qi::lexeme[(ascii::alnum >> -qi::char_('\'') >> *qi::char_("a-zA-Z_0-9") >> -(qi::char_('[') >> *ascii::digit >> qi::char_(']')))];
		vgarray %= qi::lexeme[(qi::char_('[') >> *ascii::digit >> qi::char_(':') >> *ascii::digit >> qi::char_(']'))];
		input_names %= identifier % ',';
		ios %= (qi::string("input") | qi::string("output")) >> ((identifier) % ',' | (vgarray >> identifier % ',')) >> ";";
		wire %= qi::string("wire") >> ((identifier % ',') | (vgarray >> identifier % ',')) >> ";";
		wires = *wire;	 // you may notice that = and not %= is used, that means the wires are parsed but they are not saved for later processing
		inputs_and_outputs %= *ios;
		input_pair %= qi::lit(".") >> identifier >> qi::lit("(") >> identifier >> qi::lit(")");
		gate_input_config %= qi::lit("(") >> input_pair % ',' >> qi::lit(")");
		gates %= *(identifier - "endmodule" >> identifier >> gate_input_config >> ";");

		object %= "module" >> identifier >> "(" >> input_names >> ");" >>
			inputs_and_outputs >>
			wires >>
			gates >> "endmodule";


	}
	boost::spirit::qi::rule<Iterator, circuit(), Skipper> object;
	boost::spirit::qi::rule<Iterator, std::pair<std::string, std::string>(), Skipper>input_pair;
	boost::spirit::qi::rule<Iterator, vec_type(), Skipper>gate_input_config;
	boost::spirit::qi::rule<Iterator, std::string(), Skipper>identifier, vgarray;
	boost::spirit::qi::rule<Iterator, std::vector<gate>(), Skipper> gates;
	boost::spirit::qi::rule<Iterator, std::vector<circuit_io_object>(), Skipper>inputs_and_outputs;
	boost::spirit::qi::rule<Iterator, circuit_io_object(), Skipper>ios;
	boost::spirit::qi::rule<Iterator, Skipper>wires;
	boost::spirit::qi::rule<Iterator, std::vector<std::string>(), Skipper>input_names, wire, verilogarray;
};


// function creates a calls the Boost parsing functions and parses the netlist
// into the a circuit data-structure using the above parser, filename is name of the netlist
// If parsign fails the failing message may NOT point to the correct line where the parsing failed
// but maybe the message is useful
circuit parse_netlist(const std::string &filename)
{
	circuit c;
	std::cout << "Circuit Name: " << filename << std::endl;
	std::ifstream input_file(filename.c_str());
	std::string filestream((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
	vgdata<std::string::const_iterator, boost::spirit::ascii::space_type>v;
	if (!boost::spirit::qi::phrase_parse(std::cbegin(filestream), std::cend(filestream), v, boost::spirit::ascii::space, c))
	{
		std::cout << "Parsing Failed." << std::endl;
		if (c.circuit_gates.size() > 0)
			std::cout << "Last gatename that could be parsed: " << c.circuit_gates[c.circuit_gates.size() - 1].gate_name << std::endl;
		else std::cout << "Parsing failed somewhere in the header of the file." << std::endl;

		std::cout << "Press Enter to exit Program" << std::endl;
		std::cin.get();
		//return 0;
	}
	else std::cout << "Parsing OK" << std::endl;

	return c;

}


/*
	From here the Stilefile-Parser part starts.
	The Stilefile was used to verify that the logic simulation was correct but it is not used in the final program,
	the parsing works almost in the same way the netlist does
*/
typedef std::pair<std::string,std::vector<std::string>> stilfile_signal;

struct stil_file_signal
{
	 stilfile_signal signal;
	 std::string comment_number;
};
BOOST_FUSION_ADAPT_STRUCT(stil_file_signal,(stilfile_signal, signal)
										   ( std::string, comment_number))

struct ff_order
{
	std::string circuit_name;
	std::string ff_name;
	std::string input_name;
};

BOOST_FUSION_ADAPT_STRUCT(ff_order,(std::string, circuit_name)
										   (std::string, ff_name)
										   (std::string, input_name))

typedef std::pair<std::string,std::string > input_and_pattern;
typedef std::vector<std::vector<input_and_pattern> > vec_of_vec_io;

struct allclock_capture_info
{
	std::string input_name;
	std::string input_pattern;
	std::string output_name;
	std::string output_pattern;
};

BOOST_FUSION_ADAPT_STRUCT(allclock_capture_info,(std::string, input_name)
											    (std::string, input_pattern)
												(std::string, output_name)
												(std::string, output_pattern)
						)

typedef std::pair<unsigned int,std::string > input_and_pattern_so_si;
struct pattern_information
{
	unsigned int pattern_number;
	std::vector<input_and_pattern_so_si> input_pattern_maping_test_so;
	std::vector<input_and_pattern_so_si> input_pattern_maping_test_si;
	vec_of_vec_io  multi_and_allclock_info;

};
BOOST_FUSION_ADAPT_STRUCT(pattern_information,(unsigned int, pattern_number)
											    (std::vector<input_and_pattern_so_si>, input_pattern_maping_test_so)
												(std::vector<input_and_pattern_so_si>, input_pattern_maping_test_si)
												(vec_of_vec_io, multi_and_allclock_info))


struct scan_chain
{
	unsigned int scan_chain_number;
	std::vector<ff_order> scan_structure;
};
BOOST_FUSION_ADAPT_STRUCT(scan_chain,(unsigned int, scan_chain_number)
									 (std::vector<ff_order>, scan_structure))


struct unload_data
{
	unsigned int number_of_unload_pattern;
	std::vector<input_and_pattern_so_si> test_so;
};
BOOST_FUSION_ADAPT_STRUCT(unload_data,(unsigned int, number_of_unload_pattern)
									 (std::vector<input_and_pattern_so_si>, test_so))




struct stilfile_datastruct
{
	std::vector<stil_file_signal>      signal_groups;
	std::vector<scan_chain>            scan_chains;
	std::vector<pattern_information>   pi;
	unload_data                       uload_data;
};
BOOST_FUSION_ADAPT_STRUCT(stilfile_datastruct,(std::vector<stil_file_signal>, signal_groups)
											  (std::vector<scan_chain>,scan_chains)
											  (std::vector<pattern_information>,pi)
											  (unload_data,uload_data))


template <typename Iterator, typename Skipper>
struct stilfileparser : boost::spirit::qi::grammar<Iterator, stilfile_datastruct(), Skipper>
{
	stilfileparser(): stilfileparser::base_type(object)
	{
		using namespace boost::spirit;
		identifier         %=  qi::lexeme[ qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9[]") ]; 
		binary_number      %=  *(qi::char_("0") | qi::char_("1") | qi::char_("P") | qi::char_("H") | qi::char_("L") | qi::char_("T"));
		signal_list        %=  ((qi::lit("\"") >> identifier >> qi::lit("\"")) % '+');
		signal_group_entry %=  qi::lit("\"") >> identifier >> qi::lit("\"") >> qi::lit("=") >> qi::lit("'") >> signal_list >> qi::lit("'") >> qi::lit(";"); 
		signal_groups      %=  *(signal_group_entry >> qi::lit("//") >> qi::lit("#signals=") >> *qi::char_("0-9")); 
		scan_cells         %=  *(qi::lit("\"") >> identifier >> qi::lit(".") >> identifier >> qi::lit(".") >>  identifier >> qi::lit("\""));

		key_and_value      %=  *(qi::lit("\"") >> identifier >> qi::lit("\"") >> qi::lit("=") >> binary_number >> qi::lit(";"));

		mac                %=  *(qi::lit("Call") >> qi::lit("\"") >> (qi::lit("multiclock_capture")|qi::lit("allclock_capture")|qi::lit("allclock_launch")) >> qi::lit("\"") >> qi::lit("{") >>
							   key_and_value >> qi::lit("}"));


		key_and_value_se      %=  *(qi::lit("\"") >> qi::lit("test_so") >> -boost::spirit::uint_  >> qi::lit("\"") >> qi::lit("=") >> binary_number >> qi::lit(";"));
		key_and_value_si      %=  *(qi::lit("\"") >> qi::lit("test_si") >> -boost::spirit::uint_  >> qi::lit("\"") >> qi::lit("=") >> binary_number >> qi::lit(";"));

		patterns           %= *(qi::lit("\"") >> qi::lit("pattern") >> boost::spirit::uint_ >> qi::lit("\"") >> qi::lit(":") >> 
							   qi::lit("Call") >> qi::lit("\"") >> qi::lit("load_unload") >> qi::lit("\"") >> qi::lit("{") >> key_and_value_se >> key_and_value_si >> qi::lit("}") >> mac >>
							   -qi::lit("Ann {* fast_sequential *}"));

		scanchains         %=  *(qi::lit("ScanChain") >> qi::lit("\"") >> boost::spirit::uint_  >> qi::lit("\"") >> qi::lit("{") >> qi::lit("ScanCells") >>scan_cells >> qi::lit(";"));
 
		unload_data1        %=  qi::lit("\"") >> qi::lit("end") >> -boost::spirit::uint_ >> qi::lit("unload") >> qi::lit("\"") >> qi::lit(":") >>
							   qi::lit("Call") >> qi::lit("\"") >> qi::lit("load_unload") >> qi::lit("\"") >> qi::lit("{") >> key_and_value_se >> qi::lit("}");

		object             %=  qi::lit("SignalGroups") >> qi::lit("{") >> signal_groups >> qi::lit("ScanStructures") >> qi::lit("{") >> scanchains >> patterns >> unload_data1;
	
		
	}
	boost::spirit::qi::rule<Iterator, stilfile_datastruct(),       Skipper> object;
	boost::spirit::qi::rule<Iterator, unload_data(),               Skipper> unload_data1;
	boost::spirit::qi::rule<Iterator,std::vector<pattern_information>(),Skipper>patterns;
	boost::spirit::qi::rule<Iterator, std::vector<input_and_pattern>(), Skipper> key_and_value;
	boost::spirit::qi::rule<Iterator, std::vector<input_and_pattern_so_si>(), Skipper> key_and_value_se;
	boost::spirit::qi::rule<Iterator, std::vector<input_and_pattern_so_si>(), Skipper> key_and_value_si;
	boost::spirit::qi::rule<Iterator, vec_of_vec_io(), Skipper> mac;
	boost::spirit::qi::rule<Iterator, std::vector<scan_chain>(), Skipper> scanchains; 
	boost::spirit::qi::rule<Iterator, std::vector<stil_file_signal>(),Skipper>signal_groups;
	boost::spirit::qi::rule<Iterator, std::vector<std::string>(), Skipper> signal_list;
	boost::spirit::qi::rule<Iterator, std::vector<ff_order>(), Skipper>    scan_cells;
	boost::spirit::qi::rule<Iterator, stilfile_signal()  ,Skipper> signal_group_entry;
	boost::spirit::qi::rule<Iterator, std::string(), Skipper> binary_number;
	boost::spirit::qi::rule<Iterator, std::string(), Skipper> identifier;
};

stilfile_datastruct parse_stilfile(const std::string &stilfile_filename)
{
	
	stilfile_datastruct s;
	std::ifstream stil_file(stilfile_filename.c_str());
	struct stilfile
	{
		std::string section_name;
		std::vector<std::string>section_info;
	};
	stilfile signal_groups;
	signal_groups.section_name = "SignalGroups";
	signal_groups.section_info.push_back("_pi");
	signal_groups.section_info.push_back("_po");
	/*stilfile scan_chain;
	scan_chain.section_name = "ScanStructures ";
	scan_chain.section_info.push_back("ScanCells");*/
	std::vector<stilfile>stilfilekeywords;
	stilfilekeywords.push_back(signal_groups);
	//stilfilekeywords.push_back(scan_chain);
	std::string line = "";
	std::string stilfile_data;
//	std::cout << "Reserve string size: " << 730122668 << std::endl;
	stilfile_data.reserve(730122668);
	//std::vector<std::string>stilfiledatavec;
	// Read Signal Group information
	for(auto keyword_it = stilfilekeywords.cbegin();keyword_it!=stilfilekeywords.cend();keyword_it++)
	{
		std::getline (stil_file,line);
		while(line.find((*keyword_it).section_name)==std::string::npos) std::getline (stil_file,line);
		if(line.find((*keyword_it).section_name) != std::string::npos)
		{
			stilfile_data+=std::move(line);
			//stilfiledatavec.push_back(line);
			for(auto it2 = (*keyword_it).section_info.cbegin();it2!=(*keyword_it).section_info.cend();it2++)
			{
			//	std::getline (stil_file,line);
				while(line.find(*it2)==std::string::npos)std::getline (stil_file,line);
				if(line.find(*it2) != std::string::npos)
				{
					stilfile_data+=std::move(line);
					while(line.find(";") == std::string::npos)
					{
						std::getline (stil_file,line);
						stilfile_data+=std::move(line);
					}
				}
			}
	    }
	}
	// Get the ScanChain Information
	std::getline (stil_file,line);
	while(true && stil_file.good())
	{
		std::getline (stil_file,line);
		if(line.find("ScanStructures")!=std::string::npos)
		{
			stilfile_data += std::move(line);
			while(line.find("PatternBurst")==std::string::npos)
			{
				std::getline (stil_file,line);
				if (line.find("ScanChain") != std::string::npos)stilfile_data += std::move(line);
				std::getline (stil_file,line);
				if(line.find("ScanCells")!=std::string::npos)
				{
					stilfile_data+=std::move(line);
					while(line.find(";") == std::string::npos)
					{
						std::getline (stil_file,line);
						stilfile_data+=std::move(line);
						//stilfiledatavec.push_back(line);
					}
				}
				
			}
			break;
		}
	}
	// Now read the pattern and append them to the stilfile string, not good but easy
	while(std::getline (stil_file,line))
	{
		if(line.find("pattern 0")!=std::string::npos)
		{
			stilfile_data+=std::move(line);
			//stilfiledatavec.push_back(line);
			while(std::getline (stil_file,line))
			{
			//	if(line.find("end")!=std::string::npos)break;
				 stilfile_data+=std::move(line);
			}

		}
	}

	//std::cout << "Size of Input string: " << stilfile_data.length() << std::endl;
	stilfileparser<std::string::const_iterator,boost::spirit::ascii::space_type> sfileparser;

	std::cout << "Stilfile: " << stilfile_filename << std::endl;
 	if(!boost::spirit::qi::phrase_parse(stilfile_data.cbegin(),stilfile_data.cend(),sfileparser,boost::spirit::ascii::space,s))
	{	
			std::cout << "Parsing Failed." << std::endl;
			std::cout << "Press Enter to exit Program" << std::endl;
			std::cin.get();
			//return 0;
	}
	else std::cout << "Stil file Parsing OK" << std::endl;

	return s;

}


#endif