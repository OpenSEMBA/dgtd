#include <string>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <omp.h>

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
		
	#ifndef NDEBUG
	std::cout << "PID " << getpid() << " ready to be attached. Press Enter to continue...\n";
	std::cin.get();
	#endif

	if (argc == 1) {
		std::cerr << "Missing arguments. Use argument -h to print help menu.";
		return 1;
	}
	
	std::string inputFilePath;
    std::string deviceConfig{ "cpu" };
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-i" && i + 1 < argc) {
			inputFilePath = argv[++i]; 
		}
		else if (arg == "-h") {
			printHelpArgument();
			return 2;
		}
		else if ((arg == "--device" || arg == "-d") && i + 1 < argc)
		{
			std::string devtype = argv[i+1];
			if (devtype == "cpu" || devtype == "omp" || devtype == "cuda" ){
				deviceConfig = devtype;
				++i;
			}
			else{
				throw std::runtime_error("Available device strings are \"cpu\", \"omp\" or \"cuda\"");
			}
		}
	}

    mfem::Device device(deviceConfig.c_str());
    device.Print();

	mfem::Mpi::Init(argc, argv);
	mfem::Hypre::Init();

	const std::string s_json = ".json";
	const std::string::size_type input_msh = inputFilePath.find(s_json);
	if (input_msh == std::string::npos)
	{
		std::cerr << "Input File is not a .json file." << std::endl;
		return 3;
	}

	auto solver = maxwell::driver::buildSolverJson(inputFilePath, false);

	solver.run();

	if (mfem::Mpi::WorldRank() == 0){
		std::cout << "Solver has finished performing its operations." << std::endl;
		std::cout << "Program will end now." << std::endl;
	}
	
	mfem::Mpi::Finalize();

	return 0;
}