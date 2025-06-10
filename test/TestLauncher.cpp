#include <string>
#include <gtest/gtest.h>

#include "mfem.hpp"

int main(int argc, char** argv) 
{

    std::cout << "PID " << getpid() << " ready to be attached. Press Enter to continue...\n";
    std::cin.get();

    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    std::vector<std::string> args(argc);
    for (int i = 0; i < argc; ++i) {
        args[i] = std::string(argv[i]);
        std::cout << "Argument #" << i << ": " << args[i] << std::endl;
    }

    std::string deviceConfig{ "cpu" };
    for (const auto& arg : args) {
        std::string prefix{ "--device=" };
        if (!arg.compare(0, prefix.size(), prefix)) {
            deviceConfig = arg.substr(prefix.size());
        }
    }

    mfem::Device device(deviceConfig.c_str());
    device.Print();

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    mfem::Mpi::Finalize();
    return result;
}