#include <gtest/gtest.h>
#include "components/SIBCPoleResidueFitting.cpp"

using namespace maxwell;

using Pole = Complex;
using Residue = Complex;

class SIBCFittingTest : public ::testing::Test {
};

TEST_F(SIBCFittingTest, copperCase)
{
    // Case 1: Copper (Standard Conductor / Skin Effect)
    Material mat(1.0, 1.0, 5.8e7);
    SGBCProperties props;
    props.layers.push_back({mat, 100e-6});

    auto result = findOptimalPoles(props, 1e6, 1e9, 500, 10, 0.02);

    std::vector<std::pair<Pole, Residue>> expected_data = {
        {{-1.04222e+07,  3.96072e+06}, {-1305.71,     603.538    }},
        {{-2.08739e+07, -3.35912e+07}, { 10482.4,    -33670.0    }},
        {{-4.76906e+07, -2.30975e+08}, { 31824.0,    -52752.4    }},
        {{-9.37651e+07, -9.15085e+08}, { 2.82371e+06,-3.48931e+06}},
        {{-1.91577e+08, -3.12717e+09}, {-1.30802e+07, 1.51323e+07}},
        {{-5.1975e+08,  -1.21592e+10}, { 2.8349e+08, -3.17943e+08}},
        {{-4.53152e+09, -1.28364e+11}, {-2.99485e+09, 3.3268e+09}}
    };

    auto expected_size = expected_data.size();
    ASSERT_EQ(result.poles.size(), expected_size);

    // --- COMBINE RESULT INTO PAIRS FOR SORTING ---
    std::vector<std::pair<Pole, Residue>> result_data;
    result_data.reserve(expected_size);
    for(int i = 0; i < expected_size; ++i) {
        result_data.push_back({result.poles[i], result.residues[i]});
    }

    // --- SORT BOTH LISTS BY POLE REAL PART ---
    auto sortPair = [](const std::pair<Complex, Complex>& a, const std::pair<Complex, Complex>& b) {
        return a.first.real() < b.first.real(); 
    };

    std::sort(expected_data.begin(), expected_data.end(), sortPair);
    std::sort(result_data.begin(), result_data.end(), sortPair);

// --- ROBUST COMPARISON ---
    double rel_tol = 0.05; // 5%
    double abs_tol = 1e-12; // Floor for zeros

    for (size_t i = 0; i < expected_size; ++i) {
        Pole p_res = result_data[i].first;
        Pole p_exp = expected_data[i].first;
        Residue r_res = result_data[i].second;
        Residue r_exp = expected_data[i].second;

        // Check Poles (Vector Distance)
        double p_diff = std::abs(p_res - p_exp);
        double p_mag  = std::abs(p_exp);
        EXPECT_TRUE(p_diff < (p_mag * rel_tol + abs_tol)) 
            << "Pole Mismatch at index " << i 
            << "\n  Actual: " << p_res 
            << "\n  Expect: " << p_exp;

        double r_diff = std::abs(r_res - r_exp);
        double r_mag  = std::abs(r_exp);
        EXPECT_TRUE(r_diff < (r_mag * rel_tol + abs_tol))
            << "Residue Mismatch at index " << i 
            << "\n  Actual: " << r_res 
            << "\n  Expect: " << r_exp;
    }
}

TEST_F(SIBCFittingTest, teflonCase)
{
    // Case 2: Teflon (Good Dielectric / Capacitive)
    Material mat(2.1, 1.0, 1e-12);
    SGBCProperties props;
    props.layers.push_back({mat, 100e-6});

    auto result = findOptimalPoles(props, 1e6, 1e9, 500, 10, 0.02);

    std::vector<std::pair<Pole, Residue>> expected_data = {
        {{-0.41068719, -14.54353967},  {5.37825557e+14, -6154121.42498704}}
    };

    auto expected_size = expected_data.size();
    ASSERT_EQ(result.poles.size(), expected_size);

    // --- COMBINE RESULT INTO PAIRS FOR SORTING ---
    std::vector<std::pair<Pole, Residue>> result_data;
    result_data.reserve(expected_size);
    for(int i = 0; i < expected_size; ++i) {
        result_data.push_back({result.poles[i], result.residues[i]});
    }

    // --- SORT BOTH LISTS BY POLE REAL PART ---
    auto sortPair = [](const std::pair<Pole, Residue>& a, const std::pair<Pole, Residue>& b) {
        return a.first.real() < b.first.real(); 
    };

    std::sort(expected_data.begin(), expected_data.end(), sortPair);
    std::sort(result_data.begin(), result_data.end(), sortPair);

    // --- ROBUST COMPARISON ---
    double rel_tol = 0.05; // 5%
    double abs_tol = 1e-12; // Floor for zeros

    for (size_t i = 0; i < expected_size; ++i) {
        Pole p_res = result_data[i].first;
        Pole p_exp = expected_data[i].first;
        Residue r_res = result_data[i].second;
        Residue r_exp = expected_data[i].second;

        // Check Poles (Vector Distance)
        double p_diff = std::abs(p_res - p_exp);
        double p_mag  = std::abs(p_exp);
        EXPECT_TRUE(p_diff < (p_mag * rel_tol + abs_tol)) 
            << "Pole Mismatch at index " << i 
            << "\n  Actual: " << p_res 
            << "\n  Expect: " << p_exp;

        double r_diff = std::abs(r_res - r_exp);
        double r_mag  = std::abs(r_exp);
        EXPECT_TRUE(r_diff < (r_mag * rel_tol + abs_tol))
            << "Residue Mismatch at index " << i 
            << "\n  Actual: " << r_res 
            << "\n  Expect: " << r_exp;
    }
}

TEST_F(SIBCFittingTest, saltwaterCase)
{
    // Case 2: Saltwater (Lossy Dielectric / Crossover)
    Material mat(81.0, 1.0, 4.0);
    SGBCProperties props;
    props.layers.push_back({mat, 1e-2});

    auto result = findOptimalPoles(props, 1e6, 1e9, 500, 10, 0.02);

    std::vector<std::pair<Pole, Residue>> expected_data = {
        {{-5.57089e+09, 1.96485e+07}, {1.61715e+11,-1.65994e+10}},
        {{-2.77898e+09, 1.00573e+10}, {1.10705e+11, 6.475e+10  }},
        {{-4.52649e+09,-9.82523e+09}, {2.83188e+11,-2.66864e+11}},
        {{-1.63419e+10, 2.0683e+10},  {6.91558e+11, 1.91555e+11}}
    };

    auto expected_size = expected_data.size();
    ASSERT_EQ(result.poles.size(), expected_size);

    // --- COMBINE RESULT INTO PAIRS FOR SORTING ---
    std::vector<std::pair<Pole, Residue>> result_data;
    result_data.reserve(expected_size);
    for(int i = 0; i < expected_size; ++i) {
        result_data.push_back({result.poles[i], result.residues[i]});
    }

    // --- SORT BOTH LISTS BY POLE REAL PART ---
    auto sortPair = [](const std::pair<Pole, Residue>& a, const std::pair<Pole, Residue>& b) {
        return a.first.real() < b.first.real(); 
    };

    std::sort(expected_data.begin(), expected_data.end(), sortPair);
    std::sort(result_data.begin(), result_data.end(), sortPair);

    // --- ROBUST COMPARISON ---
    double rel_tol = 0.05; // 5%
    double abs_tol = 1e-12; // Floor for zeros

    for (size_t i = 0; i < expected_size; ++i) {
        Pole p_res = result_data[i].first;
        Pole p_exp = expected_data[i].first;
        Residue r_res = result_data[i].second;
        Residue r_exp = expected_data[i].second;

        // Check Poles (Vector Distance)
        double p_diff = std::abs(p_res - p_exp);
        double p_mag  = std::abs(p_exp);
        EXPECT_TRUE(p_diff < (p_mag * rel_tol + abs_tol)) 
            << "Pole Mismatch at index " << i 
            << "\n  Actual: " << p_res 
            << "\n  Expect: " << p_exp;

        double r_diff = std::abs(r_res - r_exp);
        double r_mag  = std::abs(r_exp);
        EXPECT_TRUE(r_diff < (r_mag * rel_tol + abs_tol))
            << "Residue Mismatch at index " << i 
            << "\n  Actual: " << r_res 
            << "\n  Expect: " << r_exp;
    }
}