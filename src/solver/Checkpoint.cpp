#include "Checkpoint.h"
#include "ProbesManager.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

#include <mpi.h>

namespace maxwell {
namespace {

constexpr char kStateMagic[4] = {'D', 'G', 'T', 'C'};

uint32_t rotr(uint32_t x, uint32_t n)
{
    return (x >> n) | (x << (32u - n));
}

void sha256Block(uint32_t h[8], const uint8_t block[64])
{
    static const uint32_t k[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
        0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
        0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
        0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
        0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
        0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
        0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
    };

    uint32_t w[64];
    for (int i = 0; i < 16; ++i) {
        w[i] = (static_cast<uint32_t>(block[i * 4]) << 24)
             | (static_cast<uint32_t>(block[i * 4 + 1]) << 16)
             | (static_cast<uint32_t>(block[i * 4 + 2]) << 8)
             | static_cast<uint32_t>(block[i * 4 + 3]);
    }
    for (int i = 16; i < 64; ++i) {
        const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
    uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; ++i) {
        const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const uint32_t ch = (e & f) ^ ((~e) & g);
        const uint32_t t1 = hh + S1 + ch + k[i] + w[i];
        const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        const uint32_t t2 = S0 + maj;
        hh = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
    h[5] += f;
    h[6] += g;
    h[7] += hh;
}

class Sha256 {
public:
    Sha256()
    {
        h_[0] = 0x6a09e667u;
        h_[1] = 0xbb67ae85u;
        h_[2] = 0x3c6ef372u;
        h_[3] = 0xa54ff53au;
        h_[4] = 0x510e527fu;
        h_[5] = 0x9b05688cu;
        h_[6] = 0x1f83d9abu;
        h_[7] = 0x5be0cd19u;
    }

    void update(const void* data, std::size_t len)
    {
        const auto* bytes = static_cast<const uint8_t*>(data);
        bit_len_ += static_cast<uint64_t>(len) * 8u;
        while (len > 0) {
            const std::size_t n = std::min(len, static_cast<std::size_t>(64 - tail_));
            std::memcpy(tail_buf_ + tail_, bytes, n);
            tail_ += static_cast<int>(n);
            bytes += n;
            len -= n;
            if (tail_ == 64) {
                sha256Block(h_, tail_buf_);
                tail_ = 0;
            }
        }
    }

    std::string final()
    {
        uint8_t block[64];
        std::memcpy(block, tail_buf_, static_cast<std::size_t>(tail_));
        block[tail_] = 0x80;
        if (tail_ >= 56) {
            std::memset(block + tail_ + 1, 0, static_cast<std::size_t>(64 - (tail_ + 1)));
            sha256Block(h_, block);
            std::memset(block, 0, 56);
        } else {
            std::memset(block + tail_ + 1, 0, static_cast<std::size_t>(56 - (tail_ + 1)));
        }
        for (int i = 0; i < 8; ++i) {
            block[63 - i] = static_cast<uint8_t>(bit_len_ >> (8 * i));
        }
        sha256Block(h_, block);

        std::ostringstream out;
        out << std::hex << std::setfill('0');
        for (uint32_t word : h_) {
            out << std::setw(8) << word;
        }
        return out.str();
    }

private:
    uint32_t h_[8]{};
    uint8_t tail_buf_[64]{};
    int tail_ = 0;
    uint64_t bit_len_ = 0;
};

std::string sha256Bytes(const void* data, std::size_t len)
{
    Sha256 ctx;
    ctx.update(data, len);
    return ctx.final();
}

void checkSha256()
{
    static const bool ok = [] {
        const std::string empty = sha256Bytes(nullptr, 0);
        const char abc[] = {'a', 'b', 'c'};
        const std::string abc_hash = sha256Bytes(abc, 3);
        return empty == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
            && abc_hash == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
    }();
    if (!ok) {
        throw std::runtime_error("SHA-256 self-check failed.");
    }
}

void fsyncPath(const std::filesystem::path& path, bool directory)
{
    const int flags = directory ? (O_RDONLY | O_DIRECTORY) : O_RDONLY;
    const int fd = ::open(path.c_str(), flags);
    if (fd >= 0) {
        ::fsync(fd);
        ::close(fd);
    }
}

void durableClose(std::ofstream& out, const std::filesystem::path& path)
{
    out.flush();
    out.close();
    fsyncPath(path, false);
}

void writePod(std::ostream& out, const void* data, std::size_t bytes)
{
    out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out) {
        throw std::runtime_error("Failed while writing a checkpoint file.");
    }
}

void readPod(std::istream& in, void* data, std::size_t bytes)
{
    in.read(static_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in) {
        throw std::runtime_error("Checkpoint file is truncated.");
    }
}

int32_t toI32(int value, const char* what)
{
    if (value < 0 || static_cast<std::int64_t>(value) > static_cast<std::int64_t>(INT32_MAX)) {
        throw std::runtime_error(std::string(what) + " does not fit in the checkpoint header.");
    }
    return static_cast<int32_t>(value);
}

std::string checkpointIndexName(int index)
{
    std::ostringstream name;
    name << "checkpoint_" << std::setw(6) << std::setfill('0') << index;
    return name.str();
}

int parseCheckpointIndex(const std::string& name)
{
    constexpr const char* kPrefix = "checkpoint_";
    constexpr std::size_t kPrefixLen = 11;
    if (name.size() <= kPrefixLen || name.compare(0, kPrefixLen, kPrefix) != 0) {
        return -1;
    }
    int index = 0;
    for (std::size_t i = kPrefixLen; i < name.size(); ++i) {
        if (name[i] < '0' || name[i] > '9') {
            return -1;
        }
        index = index * 10 + (name[i] - '0');
    }
    return index;
}

void throwIfFailed(int failed, int rank, const std::string& error)
{
    int global = 0;
    MPI_Allreduce(&failed, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (!global) {
        return;
    }
    if (!error.empty()) {
        std::cerr << "Rank " << rank << ": " << error << std::endl;
    }
    if (!error.empty()) {
        throw std::runtime_error(error);
    }
    throw std::runtime_error("Checkpoint write failed on another rank.");
}

nlohmann::json cursorsToJson(const CheckpointManifest& meta)
{
    nlohmann::json exporters = nlohmann::json::array();
    for (const auto& c : meta.exporters) {
        exporters.push_back({
            {"name", c.name},
            {"save_count", c.save_count},
            {"next_save_time", c.next_save_time},
            {"dt_save", c.dt_save},
            {"initialized", c.initialized},
            {"finished", c.finished}
        });
    }
    nlohmann::json mor = nlohmann::json::array();
    for (const auto& c : meta.mor) {
        mor.push_back({
            {"name", c.name},
            {"save_count", c.save_count},
            {"next_save_time", c.next_save_time},
            {"dt_save", c.dt_save},
            {"initialized", c.initialized}
        });
    }
    return {{"exporters", exporters}, {"mor", mor}};
}

template <typename Cursor>
Cursor readCursor(const nlohmann::json& item)
{
    Cursor c;
    c.name = item.at("name").get<std::string>();
    c.save_count = item.at("save_count").get<int>();
    c.next_save_time = item.at("next_save_time").get<double>();
    c.dt_save = item.at("dt_save").get<double>();
    c.initialized = item.at("initialized").get<bool>();
    return c;
}

}  // namespace

std::string sha256File(const std::filesystem::path& path)
{
    checkSha256();
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot hash " + path.string());
    }
    Sha256 ctx;
    char buf[1 << 16];
    while (in) {
        in.read(buf, sizeof(buf));
        const auto n = in.gcount();
        if (n > 0) {
            ctx.update(buf, static_cast<std::size_t>(n));
        }
    }
    return ctx.final();
}

CheckpointSchedule checkpointSchedule(
    double time,
    double dt,
    double final_time,
    double percent,
    int next_mark)
{
    CheckpointSchedule out;
    out.next_mark = next_mark;
    if (!(percent > 0.0) || !(final_time > 0.0)) {
        return out;
    }
    double frac = time / final_time;
    if (final_time - time <= 1e-8 * std::max(dt, 1e-30)) {
        frac = 1.0;
    }
    const double step = percent / 100.0;
    if (!(step > 0.0)) {
        return out;
    }
    const int crossed = static_cast<int>(std::floor((frac + 1e-9) / step));
    if (crossed >= next_mark) {
        out.due = true;
        out.next_mark = crossed + 1;
    }
    return out;
}

bool sameCaseExceptCheckpointPercent(
    const std::filesystem::path& live,
    const std::filesystem::path& saved)
{
    try {
        std::ifstream live_file(live);
        std::ifstream saved_file(saved);
        if (!live_file || !saved_file) {
            return false;
        }
        nlohmann::json live_json;
        nlohmann::json saved_json;
        live_file >> live_json;
        saved_file >> saved_json;
        if (live_json.contains("solver_options") && live_json["solver_options"].is_object()) {
            live_json["solver_options"].erase("checkpoint_percent");
        }
        if (saved_json.contains("solver_options") && saved_json["solver_options"].is_object()) {
            saved_json["solver_options"].erase("checkpoint_percent");
        }
        return live_json == saved_json;
    } catch (const std::exception&) {
        return false;
    }
}

std::filesystem::path checkpointRoot(const std::string& caseName)
{
    return std::filesystem::path(getSimulationCaseExportPath(caseName)) / "Checkpoints";
}

namespace {

std::optional<std::filesystem::path> latestCompleteIn(const std::filesystem::path& root)
{
    if (!std::filesystem::exists(root)) {
        return std::nullopt;
    }
    int best = -1;
    std::filesystem::path best_dir;
    for (const auto& entry : std::filesystem::directory_iterator(root)) {
        if (!entry.is_directory()) {
            continue;
        }
        const int index = parseCheckpointIndex(entry.path().filename().string());
        if (index < 0 || !std::filesystem::exists(entry.path() / "COMPLETE")) {
            continue;
        }
        if (index > best) {
            best = index;
            best_dir = entry.path();
        }
    }
    if (best < 0) {
        return std::nullopt;
    }
    return best_dir;
}

}  // namespace

std::optional<std::filesystem::path> findLatestCompleteCheckpoint(const std::string& caseName)
{
    return latestCompleteIn(checkpointRoot(caseName));
}

std::optional<std::filesystem::path> findForeignCompleteCheckpoint(const std::string& caseName)
{
    const auto ours = checkpointRoot(caseName);
    const std::filesystem::path base("exports/SimulationData");
    if (!std::filesystem::exists(base)) {
        return std::nullopt;
    }
    for (const auto& mode : std::filesystem::directory_iterator(base)) {
        if (!mode.is_directory()) {
            continue;
        }
        const auto root = mode.path() / caseName / "Checkpoints";
        if (root == ours) {
            continue;
        }
        if (const auto found = latestCompleteIn(root)) {
            return found;
        }
    }
    return std::nullopt;
}

std::string foreignCheckpointMessage(const std::filesystem::path& checkpoint_dir)
{
    std::string msg = "A complete checkpoint exists at " + checkpoint_dir.string();
    try {
        const auto manifest = readManifest(checkpoint_dir / "manifest.json");
        msg += " (" + std::to_string(manifest.world_size) + " ranks)";
    } catch (const std::exception&) {
    }
    msg += ". This job is " + getRunModeTag()
        + ". Resume with the same MPI rank count and device, from the same working directory.";
    return msg;
}

CheckpointManifest readManifest(const std::filesystem::path& manifest_json)
{
    std::ifstream in(manifest_json);
    if (!in) {
        throw std::runtime_error("Cannot read checkpoint manifest " + manifest_json.string());
    }
    nlohmann::json j;
    in >> j;
    CheckpointManifest m;
    m.format_version = j.at("format_version").get<int>();
    m.json_sha256 = j.at("json_sha256").get<std::string>();
    m.mesh_sha256 = j.at("mesh_sha256").get<std::string>();
    m.world_size = j.at("world_size").get<int>();
    m.time = j.at("time").get<double>();
    m.dt = j.at("dt").get<double>();
    m.final_time = j.at("final_time").get<double>();
    m.order = j.at("order").get<int>();
    m.cycle = j.at("cycle").get<int>();
    m.next_checkpoint_mark = j.at("next_checkpoint_mark").get<int>();
    m.state_sizes = j.at("state_sizes").get<std::vector<int>>();
    m.sgbc_counts = j.at("sgbc_counts").get<std::vector<int>>();
    for (const auto& item : j.at("exporters")) {
        auto c = readCursor<ExporterCursor>(item);
        c.finished = item.at("finished").get<bool>();
        m.exporters.push_back(std::move(c));
    }
    for (const auto& item : j.at("mor")) {
        m.mor.push_back(readCursor<MorCursor>(item));
    }
    if (j.contains("elapsed_run_seconds")) {
        m.elapsed_run_seconds = j.at("elapsed_run_seconds").get<double>();
    }
    if (j.contains("state_sha256")) {
        m.state_sha256 = j.at("state_sha256").get<std::vector<std::string>>();
    }
    if (j.contains("temporal_mem_baseline")) {
        m.temporal_mem_baseline = j.at("temporal_mem_baseline").get<std::vector<std::uint64_t>>();
    }
    if (j.contains("temporal_mem_peak")) {
        m.temporal_mem_peak = j.at("temporal_mem_peak").get<std::vector<std::uint64_t>>();
    }
    if (j.contains("temporal_mem_sum")) {
        m.temporal_mem_sum = j.at("temporal_mem_sum").get<std::vector<double>>();
    }
    if (j.contains("temporal_mem_count")) {
        m.temporal_mem_count = j.at("temporal_mem_count").get<std::vector<std::uint64_t>>();
    }
    return m;
}

std::vector<int> readPartition(const std::filesystem::path& partition_bin)
{
    std::ifstream in(partition_bin, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot read checkpoint partition " + partition_bin.string());
    }
    int32_t n = 0;
    readPod(in, &n, sizeof(n));
    if (n < 0) {
        throw std::runtime_error("Checkpoint partition count is negative.");
    }
    std::vector<int32_t> raw(static_cast<std::size_t>(n));
    if (n > 0) {
        readPod(in, raw.data(), static_cast<std::size_t>(n) * sizeof(int32_t));
    }
    return std::vector<int>(raw.begin(), raw.end());
}

void readRankState(
    const std::filesystem::path& state_file,
    double* state,
    int state_size,
    std::vector<SgbcRecord>& sgbc)
{
    std::ifstream in(state_file, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot read " + state_file.string());
    }
    char magic[4];
    readPod(in, magic, 4);
    if (std::memcmp(magic, kStateMagic, 4) != 0) {
        throw std::runtime_error("Unrecognized checkpoint state file " + state_file.string());
    }
    int32_t version = 0;
    readPod(in, &version, sizeof(version));
    if (version != kCheckpointFormatVersion) {
        throw std::runtime_error(
            "Checkpoint state version " + std::to_string(version)
            + " is not supported.");
    }
    int32_t nstate = 0;
    readPod(in, &nstate, sizeof(nstate));
    if (nstate != state_size) {
        throw std::runtime_error(
            "Checkpoint state length " + std::to_string(nstate)
            + " does not match this rank's state length " + std::to_string(state_size) + ".");
    }
    if (nstate > 0) {
        readPod(in, state, static_cast<std::size_t>(nstate) * sizeof(double));
    }
    int32_t nsgbc = 0;
    readPod(in, &nsgbc, sizeof(nsgbc));
    if (nsgbc < 0) {
        throw std::runtime_error("Checkpoint SGBC count is negative.");
    }
    sgbc.clear();
    sgbc.reserve(static_cast<std::size_t>(nsgbc));
    for (int32_t i = 0; i < nsgbc; ++i) {
        SgbcRecord rec;
        int32_t tag = 0, a = 0, b = 0, nvals = 0;
        readPod(in, &tag, sizeof(tag));
        readPod(in, &a, sizeof(a));
        readPod(in, &b, sizeof(b));
        readPod(in, &nvals, sizeof(nvals));
        if (nvals < 0) {
            throw std::runtime_error("Checkpoint SGBC state length is negative.");
        }
        rec.tag = tag;
        rec.node_a = a;
        rec.node_b = b;
        rec.values.resize(static_cast<std::size_t>(nvals));
        if (nvals > 0) {
            readPod(in, rec.values.data(), static_cast<std::size_t>(nvals) * sizeof(double));
        }
        sgbc.push_back(std::move(rec));
    }
}

void writeCheckpointCollective(
    const std::string& case_name,
    const std::string& json_path,
    const std::string& mesh_path,
    const CheckpointManifest& meta,
    const int* partition,
    int partition_size,
    const double* state,
    int state_size,
    const std::vector<SgbcRecord>& sgbc)
{
    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const auto root = checkpointRoot(case_name);
    const auto partial = root / "partial";
    std::filesystem::path dest;
    std::string error;
    int failed = 0;

    if (rank == 0) {
        try {
            std::filesystem::create_directories(root);
            if (std::filesystem::exists(partial)) {
                std::filesystem::remove_all(partial);
            }
            std::filesystem::create_directories(partial);
        } catch (const std::exception& ex) {
            failed = 1;
            error = ex.what();
        }
    }
    MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (failed) {
        if (rank == 0 && !error.empty()) {
            std::cerr << error << std::endl;
        }
        throw std::runtime_error(error.empty()
            ? std::string("Checkpoint directory could not be created.")
            : error);
    }

    const auto state_path = partial / ("state.rank" + std::to_string(rank) + ".bin");
    std::string local_hash;
    try {
        std::ofstream out(state_path, std::ios::binary | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("Cannot write " + state_path.string());
        }
        Sha256 hasher;
        const auto emit = [&](const void* data, std::size_t bytes) {
            hasher.update(data, bytes);
            writePod(out, data, bytes);
        };
        emit(kStateMagic, 4);
        const int32_t version = kCheckpointFormatVersion;
        const int32_t nstate = toI32(state_size, "State length");
        const int32_t nsgbc = toI32(static_cast<int>(sgbc.size()), "SGBC count");
        emit(&version, sizeof(version));
        emit(&nstate, sizeof(nstate));
        if (nstate > 0) {
            emit(state, static_cast<std::size_t>(nstate) * sizeof(double));
        }
        emit(&nsgbc, sizeof(nsgbc));
        for (const auto& rec : sgbc) {
            const int32_t tag = toI32(rec.tag, "SGBC tag");
            const int32_t a = rec.node_a;
            const int32_t b = rec.node_b;
            const int32_t nvals = toI32(static_cast<int>(rec.values.size()), "SGBC state length");
            emit(&tag, sizeof(tag));
            emit(&a, sizeof(a));
            emit(&b, sizeof(b));
            emit(&nvals, sizeof(nvals));
            if (nvals > 0) {
                emit(rec.values.data(), static_cast<std::size_t>(nvals) * sizeof(double));
            }
        }
        local_hash = hasher.final();
        durableClose(out, state_path);
    } catch (const std::exception& ex) {
        failed = 1;
        error = ex.what();
    }
    throwIfFailed(failed, rank, error);

    std::vector<char> hash_gather(static_cast<std::size_t>(nprocs) * 65);
    char hash_local[65] = {};
    std::memcpy(hash_local, local_hash.data(), std::min(local_hash.size(), sizeof(hash_local) - 1));
    MPI_Gather(hash_local, 65, MPI_CHAR, hash_gather.data(), 65, MPI_CHAR, 0, MPI_COMM_WORLD);
    std::vector<std::string> state_sha256(static_cast<std::size_t>(nprocs));
    if (rank == 0) {
        for (int r = 0; r < nprocs; ++r) {
            state_sha256[static_cast<std::size_t>(r)] = std::string(hash_gather.data() + static_cast<std::size_t>(r) * 65);
        }
    }

    std::vector<int> state_sizes(static_cast<std::size_t>(nprocs));
    std::vector<int> sgbc_counts(static_cast<std::size_t>(nprocs));
    const int local_sgbc = static_cast<int>(sgbc.size());
    MPI_Gather(&state_size, 1, MPI_INT, state_sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_sgbc, 1, MPI_INT, sgbc_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        try {
            CheckpointManifest written = meta;
            written.format_version = kCheckpointFormatVersion;
            written.world_size = nprocs;
            written.state_sizes = std::move(state_sizes);
            written.sgbc_counts = std::move(sgbc_counts);
            written.state_sha256 = std::move(state_sha256);
            if (written.json_sha256.empty()) {
                written.json_sha256 = sha256File(json_path);
            }
            if (written.mesh_sha256.empty()) {
                written.mesh_sha256 = sha256File(mesh_path);
            }

            nlohmann::json j = {
                {"format_version", written.format_version},
                {"json_sha256", written.json_sha256},
                {"mesh_sha256", written.mesh_sha256},
                {"world_size", written.world_size},
                {"time", written.time},
                {"dt", written.dt},
                {"final_time", written.final_time},
                {"order", written.order},
                {"cycle", written.cycle},
                {"next_checkpoint_mark", written.next_checkpoint_mark},
                {"elapsed_run_seconds", written.elapsed_run_seconds},
                {"state_sizes", written.state_sizes},
                {"sgbc_counts", written.sgbc_counts},
                {"state_sha256", written.state_sha256},
                {"temporal_mem_baseline", written.temporal_mem_baseline},
                {"temporal_mem_peak", written.temporal_mem_peak},
                {"temporal_mem_sum", written.temporal_mem_sum},
                {"temporal_mem_count", written.temporal_mem_count}
            };
            const auto cursors = cursorsToJson(written);
            j["exporters"] = cursors["exporters"];
            j["mor"] = cursors["mor"];

            const auto manifest_path = partial / "manifest.json";
            {
                std::ofstream out(manifest_path);
                if (!out) {
                    throw std::runtime_error("Cannot write " + manifest_path.string());
                }
                out << j.dump(2) << '\n';
                durableClose(out, manifest_path);
            }

            const auto part_path = partial / "partition.bin";
            {
                if (partition_size > 0 && partition == nullptr) {
                    throw std::runtime_error("Checkpoint partition is missing.");
                }
                std::ofstream out(part_path, std::ios::binary | std::ios::trunc);
                if (!out) {
                    throw std::runtime_error("Cannot write " + part_path.string());
                }
                const int32_t n = toI32(partition_size, "Partition length");
                writePod(out, &n, sizeof(n));
                for (int i = 0; i < partition_size; ++i) {
                    const int32_t rank_id = toI32(partition[i], "Partition rank");
                    writePod(out, &rank_id, sizeof(rank_id));
                }
                durableClose(out, part_path);
            }

            std::filesystem::copy_file(
                json_path,
                partial / "input.json",
                std::filesystem::copy_options::overwrite_existing);
            fsyncPath(partial / "input.json", false);

            int next_index = 0;
            if (std::filesystem::exists(root)) {
                for (const auto& entry : std::filesystem::directory_iterator(root)) {
                    if (!entry.is_directory()) {
                        continue;
                    }
                    const int index = parseCheckpointIndex(entry.path().filename().string());
                    if (index > next_index) {
                        next_index = index;
                    }
                }
            }
            ++next_index;
            dest = root / checkpointIndexName(next_index);
            if (std::filesystem::exists(dest)) {
                throw std::runtime_error("Checkpoint directory already exists: " + dest.string());
            }
            std::filesystem::rename(partial, dest);
            {
                std::ofstream complete(dest / "COMPLETE");
                if (!complete) {
                    throw std::runtime_error("Cannot write COMPLETE marker in " + dest.string());
                }
                complete << "ok\n";
                durableClose(complete, dest / "COMPLETE");
            }
            fsyncPath(dest, true);

            for (const auto& entry : std::filesystem::directory_iterator(root)) {
                if (!entry.is_directory() || entry.path() == dest) {
                    continue;
                }
                const int index = parseCheckpointIndex(entry.path().filename().string());
                if (index >= 0) {
                    std::filesystem::remove_all(entry.path());
                }
            }
            fsyncPath(root, true);

            std::cout << "Checkpoint " << dest.string() << " at t = " << written.time;
            if (written.final_time > 0.0) {
                std::cout << " (" << (100.0 * written.time / written.final_time) << "%)";
            }
            std::cout << std::endl;
        } catch (const std::exception& ex) {
            failed = 1;
            error = ex.what();
        }
    }

    MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (failed) {
        if (rank == 0 && !error.empty()) {
            std::cerr << error << std::endl;
        }
        throw std::runtime_error(error.empty()
            ? std::string("Checkpoint commit failed.")
            : error);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

}  // namespace maxwell
