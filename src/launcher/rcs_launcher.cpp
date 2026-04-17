#include <iostream>
#include <string>

#include <mfem.hpp>

#include "driver/rcs_driver.h"

void printHelpArgument()
{
	std::cout << "                     OpenSEMBA/dgtd RCS Post-Processor                 "  << std::endl;
	std::cout << "___________________________________________________________________________"  << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "Command line arguments:"                                                      << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "___________________________________________________________________________"  << std::endl;
	std::cout <<																				   std::endl;
	std::cout << "-h              / Brings up this help menu."								    << std::endl;
	std::cout << "-i path/to/file / Specifies RCS input JSON file."							    << std::endl;
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
		if (arg == "-i" && i + 1 < argc) {
			inputFilePath = argv[++i];
		}
		else if (arg == "-h") {
			printHelpArgument();
			return 0;
		}
	}

	if (inputFilePath.empty()) {
		std::cerr << "No input file specified. Use -i path/to/file.json" << std::endl;
		return 1;
	}

	mfem::Mpi::Init(argc, argv);
	mfem::Hypre::Init();

	maxwell::driver::runRCSPostProcessing(inputFilePath);

	if (mfem::Mpi::WorldRank() == 0) {
		std::cout << "RCS post-processing complete." << std::endl;
	}

	mfem::Mpi::Finalize();

	return 0;
}
