# Description
This program can be used to simulate path delay fault in digital circuits as
described in https://ieeexplore.ieee.org/document/6979116

The flow of the program is as follows:
* Read a netlist (circuit data) from a file (Boost Spirit is used for parsing the files)
* Create a graph based on the data read in the previous step (the Boost Graph library is used for this part)
* Run the simulation, the simulation function is written in OpenCL and can be executed on the CPU and GPU

# Requirements
In order to compile the project, the compiler needs to support C++11. This project was developed with Visual Studio 2013 Desktop edition. But the project may also compile with the g++ or other compilers. Only standard C++ is used, no operating system specific functions are used.
 
The project makes usage of the Boost libraries and the OpenCL SDK. These two libraries need to be installed in order to compile the project.

| Library Name | Version                                      |
| ---          | ---                                          |
| Boost        | 1.55 or newer (older versions may also work) |
| OpenCL SDK   | Any version is ok                            |

	
# Short description for each file

The project consists of nine files. One “.cpp” file, seven “.h” files and one “.cl” file.

|File Name	|Short description|
| ---          | ---                                          |
|main.cpp	|The entry point of the program. Inside the main function several functions from the header files are called in order so parse the graph and set-up the simulation.|
|file_parsers.h	|In here the parsing of the netlist takes place. This part of the program makes usage of the Boost.Spirit library.|
|graph_creator.h	|In here the construction of the graph from the netlist takes place. |
|graph_structure.h	|This file contains the definition of the graph, that means how the data is stored in the graph data structure. |
|hash_maps.h	|This file defines several hash maps which are used during for creating the graph.|
|helper_functions.h	|This file defines several small functions which make it easy to work with the program, e.g it contains a function to levelize a circuit. |
|logic_simulation.h	|Contains two functions for performing a simple logic simulation on the CPU in a sequential manner.|
|pathclassification.h	|Contains the algorithm for the path classification. Note, the algorithm already contains a logic simulation. The logic simulations defined in the file logic_simualtion.h were developed for the stile file checking.|
|PathDelayFault.cl	|This file contains the OpenCL program which performs the simulation in parallel on the GPU or CPU.|

# Description of configuration variables
With the help of these variables different versions of the program can be compiled e.g. Program performs just a logic simulation, Program performs path classification on the CPU etc.

|Variable Name |Description                                   |
| ---          | ---                                          |
|CPU_LOGIC_AND_PATH_SIM	| If this variable is set to “true”, the program performs a sequential path classification on the CPU. |
|OPENCL_SIMULATION	| If this variable is set to “true” the program performs a parallel path classification on the CPU.    |
|OPENCL_GPU_PATH_CLASSIFICATION	| If this variable is set to “true” the program performs a parallel path classification on the GPU.| 
|STILEFILE_CHECK	|Read the Stilfile which is passed as a command line argument to the program and verifies the circuit. |
|LOGGING	        |If set to true the final data which is generated by the sequential path classification is logged and written into a “.dot” file. From this file an image of the graph can be generated. |
|OPENCL_LOGGING	| If set to true the final data which is generated by the parallel OpenCL version is logged and written into a files. Note, for every thread a separate “.dot” file will be generated. From each “.dot” file a graph can be generated. |


# How to run the program
In order to run the program, three parameters must be passed to the program:

```
SimpleParser.exe   netlist.vg  dummy.stil  NumberOfThreads
```

Actually only the first and last parameter are used but for the parameter in the middle a file should be passed which exists in the directory. The file in the middle must be a valid stilefile when the STILEFILE_CHECK is set to true. When the sequential version of the algorithm is used NumberOfThreads corresponds to the number of pattern which are simulated. If the any of the OpenCL versions is used the number of patterns which are simulated is: NumberOfThreads * 32 
