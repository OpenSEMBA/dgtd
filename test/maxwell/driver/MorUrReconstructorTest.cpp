#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "driver/MorUrReconstructor.h"

using namespace maxwell::driver;

namespace {

std::filesystem::path makeTempDir()
{
    const auto dir =
        std::filesystem::temp_directory_path() /
        ("mor_ur_test_" + std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    return dir;
}

} // namespace

TEST(MorUrReconstructor, urMatVecMatchesReference)
{
    constexpr std::int64_t N = 5;
    constexpr std::int64_t r = 3;
    // Column-major Ur:
    // cols: [1,2,3,4,5]^T, [0,1,0,1,0]^T, [2,2,2,2,2]^T
    std::vector<double> Ur = {
        1, 2, 3, 4, 5,
        0, 1, 0, 1, 0,
        2, 2, 2, 2, 2
    };
    mfem::Vector xr(static_cast<int>(r));
    xr[0] = 1.0;
    xr[1] = -2.0;
    xr[2] = 0.5;

    mfem::Vector x(static_cast<int>(N));
    urMatVec(Ur, N, r, xr, x);

    // x = col0 - 2*col1 + 0.5*col2
    EXPECT_DOUBLE_EQ(x[0], 1.0 - 0.0 + 1.0);
    EXPECT_DOUBLE_EQ(x[1], 2.0 - 2.0 + 1.0);
    EXPECT_DOUBLE_EQ(x[2], 3.0 - 0.0 + 1.0);
    EXPECT_DOUBLE_EQ(x[3], 4.0 - 2.0 + 1.0);
    EXPECT_DOUBLE_EQ(x[4], 5.0 - 0.0 + 1.0);
}

TEST(MorUrReconstructor, loadUrBinaryAndMetaRoundTrip)
{
    const auto dir = makeTempDir();
    const auto ur_path = dir / "Ur.bin";
    const auto meta_path = dir / "meta.json";

    constexpr std::int64_t N = 4;
    constexpr std::int64_t r = 2;
    std::vector<double> Ur = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    {
        std::ofstream ofs(ur_path, std::ios::binary);
        ofs.write(reinterpret_cast<const char*>(Ur.data()),
                  static_cast<std::streamsize>(Ur.size() * sizeof(double)));
    }
    {
        std::ofstream ofs(meta_path);
        ofs << R"({
  "N": 4,
  "r": 2,
  "dtype": "float64",
  "layout": "colmajor",
  "ur_file": "Ur.bin",
  "mor_case": "unit_test"
})";
    }

    const auto meta = loadMorUrMeta(meta_path.string());
    EXPECT_EQ(meta.N, 4);
    EXPECT_EQ(meta.r, 2);
    EXPECT_EQ(meta.layout, "colmajor");

    const auto loaded = loadUrBinary(ur_path.string(), meta.N, meta.r, meta.layout, meta.dtype);
    ASSERT_EQ(loaded.size(), Ur.size());
    for (std::size_t i = 0; i < Ur.size(); ++i) {
        EXPECT_DOUBLE_EQ(loaded[i], Ur[i]);
    }

    mfem::Vector xr(2);
    xr = 1.0;
    mfem::Vector x(4);
    urMatVec(loaded, meta.N, meta.r, xr, x);
    EXPECT_DOUBLE_EQ(x[0], 6.0);
    EXPECT_DOUBLE_EQ(x[1], 8.0);
    EXPECT_DOUBLE_EQ(x[2], 10.0);
    EXPECT_DOUBLE_EQ(x[3], 12.0);

    std::filesystem::remove_all(dir);
}

TEST(MorUrReconstructor, collectXrAcceptsTxtSuffix)
{
    const auto dir = makeTempDir();
    {
        std::ofstream ofs(dir / "xr_0.txt");
        ofs << "0.0\n2\n1.0\n2.0\n";
    }
    {
        std::ofstream ofs(dir / "xr_2.txt");
        ofs << "0.2\n2\n3.0\n4.0\n";
    }
    {
        std::ofstream ofs(dir / "xr_1"); // bare name also OK
        ofs << "0.1\n2\n5.0\n6.0\n";
    }

    const auto files = collectIndexedSnapshots(dir.string(), "xr");
    ASSERT_EQ(files.size(), 3u);
    EXPECT_EQ(files[0].first, 0);
    EXPECT_EQ(files[1].first, 1);
    EXPECT_EQ(files[2].first, 2);
    EXPECT_NE(files[0].second.find("xr_0.txt"), std::string::npos);
    EXPECT_NE(files[1].second.find("xr_1"), std::string::npos);
    EXPECT_NE(files[2].second.find("xr_2.txt"), std::string::npos);

    std::filesystem::remove_all(dir);
}

TEST(MorUrReconstructor, loadAsciiXrSnapshot)
{
    const auto dir = makeTempDir();
    const auto xr_path = dir / "xr_0.txt";
    {
        std::ofstream ofs(xr_path);
        ofs << "1.25\n3\n0.5\n-1.0\n2.0\n";
    }

    const auto [time, xr] = loadAsciiSnapshot(xr_path.string(), 3);
    EXPECT_DOUBLE_EQ(time, 1.25);
    ASSERT_EQ(xr.Size(), 3);
    EXPECT_DOUBLE_EQ(xr[0], 0.5);
    EXPECT_DOUBLE_EQ(xr[1], -1.0);
    EXPECT_DOUBLE_EQ(xr[2], 2.0);

    std::filesystem::remove_all(dir);
}
