// #include "L2ErrorAnalysis.h"
// #include <gtest/gtest.h>
// #include <filesystem>
// #include <iostream>
// #include <sstream>

// class AnalyticalStudyTests : public ::testing::Test {
// protected:
//     void SetUp() override {}
//     void TearDown() override {}
// };

// // Generic batch runner that takes a prefix and P-range
// static void runBatchSweep(const std::string& case_prefix, int p_min, int p_max, const std::string& description) 
// {
//     // Adjust this base path if your exports are stored elsewhere
//     const std::string base_export_path = "./Exports/single-core/"; 
    
//     int processed = 0;
//     int missing = 0;

//     std::cout << "\n=== Starting Batch: " << description << " ===" << std::endl;

//     for (int p = p_min; p <= p_max; ++p) {
        
//         std::stringstream ss;
//         // Construct folder name: e.g. "2D_BesselJ6_Coarse_G2_H0_P7"
//         // The prefix should include everything up to "_P"
//         ss << case_prefix << "_P" << p;
        
//         std::string case_name = ss.str();
//         std::filesystem::path full_path = std::filesystem::path(base_export_path) / case_name;

//         if (std::filesystem::exists(full_path) && std::filesystem::is_directory(full_path)) {
//             std::cout << "Processing: " << case_name << " ... ";
//             try {
//                 maxwell::L2ErrorAnalysis analysis(full_path.string());
//                 std::cout << "[DONE]" << std::endl;
//                 processed++;
//             }
//             catch (const std::exception& e) {
//                 std::cerr << "\n[ERROR] in " << case_name << ": " << e.what() << std::endl;
//             }
//         } else {
//             // Optional: Print missed files to debug paths
//             // std::cout << "Skipping (Not Found): " << case_name << std::endl;
//             missing++;
//         }
//     }

//     std::cout << "=== Batch Complete (" << description << ") ===\n";
//     std::cout << "Processed: " << processed << ", Skipped (Missing): " << missing << "\n" << std::endl;
// }

// TEST_F(AnalyticalStudyTests, Batch_BesselJ6_G1_H0_Linear_Sweeps)
// {
//     runBatchSweep("2D_BesselJ6_G1_H0", 1, 5, "Bessel J6 - Linear");
// }

// TEST_F(AnalyticalStudyTests, Batch_BesselJ6_G2_H0_Curved_Sweeps)
// {
//     runBatchSweep("2D_BesselJ6_G2_H0", 1, 5, "Bessel J6 - Curved");
// }

// // TEST_F(AnalyticalStudyTests, Batch_TM55_Resonant_Box)
// // {
// //    // This would need the old double-loop logic or manual calls
// //    // runBatchSweep("2D_Resonant_Box_TM55_H1", 1, 10, "TM55 H1");
// //    // runBatchSweep("2D_Resonant_Box_TM55_H2", 1, 10, "TM55 H2");
// //    // runBatchSweep("2D_Resonant_Box_TM55_H3", 1, 10, "TM55 H3");
// // }