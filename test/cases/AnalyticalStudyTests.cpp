#include "L2ErrorAnalysis.h"

#include <filesystem>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>

class AnalyticalStudyTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

static void runBatchSweep(const std::string& case_prefix,
                          const std::string& case_suffix,
                          int p_min,
                          int p_max,
                          const std::string& description)
{
    const std::string base_export_path = "./Exports/single-core/";

    int processed = 0;
    int missing = 0;

    std::cout << "\n=== Starting Batch: " << description << " ===" << std::endl;

    for (int p = p_min; p <= p_max; ++p) {
        std::stringstream ss;
        ss << case_prefix << "_P" << p << case_suffix;

        const std::string case_name = ss.str();
        const std::filesystem::path full_path = std::filesystem::path(base_export_path) / case_name;

        if (std::filesystem::exists(full_path) && std::filesystem::is_directory(full_path)) {
            std::cout << "Processing: " << case_name << " ... ";
            try {
                maxwell::L2ErrorAnalysis analysis(full_path.string());
                std::cout << "[DONE]" << std::endl;
                processed++;
            } catch (const std::exception& e) {
                std::cerr << "\n[ERROR] in " << case_name << ": " << e.what() << std::endl;
            }
        } else {
            missing++;
        }
    }

    std::cout << "=== Batch Complete (" << description << ") ===\n";
    std::cout << "Processed: " << processed << ", Skipped (Missing): " << missing << "\n" << std::endl;

    EXPECT_GT(processed, 0) << "No matching export folders were found for batch: " << description;
}

static void runResonantFamily(const std::string& suffix, const std::string& description)
{
    for (int h = 1; h <= 3; ++h) {
        std::stringstream prefix;
        prefix << "2D_Resonant_Box_TM55_H" << h;
        std::stringstream batch_desc;
        batch_desc << description << " H" << h;
        runBatchSweep(prefix.str(), suffix, 1, 10, batch_desc.str());
    }
}

TEST_F(AnalyticalStudyTests, Batch_TM55_Resonant_Box_Global_Closed)
{
    runResonantFamily("", "TM55 Global Closed Basis");
}

TEST_F(AnalyticalStudyTests, Batch_TM55_Resonant_Box_Hesthaven_Closed)
{
    runResonantFamily("_hesthaven", "TM55 Hesthaven Closed Basis");
}

TEST_F(AnalyticalStudyTests, Batch_TM55_Resonant_Box_Global_Open)
{
    runResonantFamily("_bleg", "TM55 Global Open Basis");
}
