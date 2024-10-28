#include <string>
#include <vector>
#include <iostream>

#include "driver/driver.h"

void printHelpArgument() 
{
	std::cout << "                               OpenSEMBA/dgtd                              "  << std::endl;
	std::cout << "___________________________________________________________________________"  << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "Command line arguments:"                                                      << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "___________________________________________________________________________"  << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "-h              / Brings up this help menu."								    << std::endl;
	std::cout << "-i path/to/file / Specifies input file to run with executable."			    << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "___________________________________________________________________________"  << std::endl;
}

int main(int argc, char** argv)
{
	
	if (argc == 1) {
		std::cerr << "Missing arguments. Use argument -h to print help menu.";
		return 1;
	}
	
	std::string inputFilePath;

	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (std::string(argv[i]) == "-i" && i + 1 < argc) {
			inputFilePath = argv[++i]; 
		}
		else if (std::string(argv[i]) == "-h") {
			printHelpArgument();
			system("pause");
			return 2;
		}
	}

	auto solver = maxwell::driver::buildSolverJson(inputFilePath, false);
	solver.run();

	std::cout << "Solver has finished performing its operations." << std::endl;
	std::cout << "Program will end now." << std::endl;

	return 0;
}