#include <string>
#include <vector>
#include <iostream>
#include <gtest/gtest.h>

#include "mfem.hpp"

int main(int argc, char** argv) 
{
    #ifndef NDEBUG
    std::cout << "PID " << getpid() << " ready to be attached. Press Enter to continue...\n";
    std::cin.get();
    #endif

    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    std::string deviceConfig{ "cpu" };
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "--device" || arg == "-d") && i + 1 < argc) {
            std::string devtype = argv[i + 1];
            if (devtype == "cpu" || devtype == "omp" || devtype == "cuda") {
                deviceConfig = devtype;
                ++i;
            } else {
                throw std::runtime_error(
                    "Available device strings are \"cpu\", \"omp\" or \"cuda\"");
            }
        }
    }

    mfem::Device device(deviceConfig.c_str());
    device.Print();

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    mfem::Mpi::Finalize();
    return result;
}
