# Description
This program can be used to simulate path delay fault simulations in digital circuits.
The program makes usage of several Boost libraries (only header libraries).

The flow of the program is as follows:
* Read a netlist (circuit data) from a file (Boost Spirit is used for parsing the files)
* Create a graph based on the data read in the previous step (the Boost Graph library is used for this part)
* Run the simulation, the simulation function is written in OpenCL and can be executed on the CPU and GPU
