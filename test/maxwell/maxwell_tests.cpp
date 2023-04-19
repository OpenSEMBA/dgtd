#include <string>
#include "gtest/gtest.h"

#include "mfem.hpp"

int main(int argc, char** argv) 
{
	std::string command_line_arg(argc == 2 ? argv[1] : "");
	
	std::vector<std::string> args(argc);
	for (auto i{ 0 }; i < argc; i++) {
		args[i] = std::string(argv[i]);
		std::cout << "Argument #" << i << ": " << args[i] << std::endl;
	}

	std::string deviceConfig{ "cpu" };
	for (const auto& arg : args) {
		std::string prefix{"--device="};
		if (!arg.compare(0, prefix.size(), prefix)) {
			deviceConfig = arg.substr(prefix.size());
		}
	}
	
	mfem::Device device(deviceConfig.c_str());
	device.Print();

	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}