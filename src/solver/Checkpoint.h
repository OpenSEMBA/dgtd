#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace maxwell {

inline constexpr int kCheckpointFormatVersion = 1;

struct ExporterCursor {
    std::string name;
    int save_count = 0;
    double next_save_time = 0.0;
    double dt_save = 0.0;
    bool initialized = false;
    bool finished = false;
};

struct MorCursor {
    std::string name;
    int save_count = 0;
    double next_save_time = 0.0;
    double dt_save = 0.0;
    bool initialized = false;
};

struct CheckpointManifest {
    int format_version = kCheckpointFormatVersion;
    std::string json_sha256;
    std::string mesh_sha256;
    int world_size = 0;
    double time = 0.0;
    double dt = 0.0;
    double final_time = 0.0;
    int order = 0;
    int cycle = 0;
    int next_checkpoint_mark = 1;
    /// Solver::run seconds already stored in earlier checkpoints. The minutes
    /// after this save and before a crash are not included.
    double elapsed_run_seconds = 0.0;
    std::vector<int> state_sizes;
    std::vector<int> sgbc_counts;
    std::vector<std::string> state_sha256;
    std::vector<std::uint64_t> temporal_mem_baseline;
    std::vector<std::uint64_t> temporal_mem_peak;
    std::vector<double> temporal_mem_sum;
    std::vector<std::uint64_t> temporal_mem_count;
    std::vector<ExporterCursor> exporters;
    std::vector<MorCursor> mor;
};

struct SgbcRecord {
    int tag = 0;
    int node_a = 0;
    int node_b = -1;
    std::vector<double> values;
};

std::string sha256File(const std::filesystem::path& path);

std::filesystem::path checkpointRoot(const std::string& caseName);

std::optional<std::filesystem::path> findLatestCompleteCheckpoint(const std::string& caseName);

/// A complete checkpoint for this case written under a different run-mode
/// directory (another -np or device). Empty when the only saves belong to
/// this job's directory.
std::optional<std::filesystem::path> findForeignCompleteCheckpoint(const std::string& caseName);

std::string foreignCheckpointMessage(const std::filesystem::path& checkpoint_dir);

CheckpointManifest readManifest(const std::filesystem::path& manifest_json);

std::vector<int> readPartition(const std::filesystem::path& partition_bin);

/// Collective over MPI_COMM_WORLD. Rank 0 writes the manifest, partition, and
/// JSON copy. Every rank writes state.rank<r>.bin. The previous complete
/// checkpoint is removed only after the new COMPLETE marker is in place.
void writeCheckpointCollective(
    const std::string& case_name,
    const std::string& json_path,
    const std::string& mesh_path,
    const CheckpointManifest& meta,
    const int* partition,
    int partition_size,
    const double* state,
    int state_size,
    const std::vector<SgbcRecord>& sgbc);

void readRankState(
    const std::filesystem::path& state_file,
    double* state,
    int state_size,
    std::vector<SgbcRecord>& sgbc);

}
