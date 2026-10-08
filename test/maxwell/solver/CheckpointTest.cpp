#include <gtest/gtest.h>

#include "solver/Checkpoint.h"
#include "solver/ProbesManager.h"

#include <fstream>
#include <mpi.h>

using namespace maxwell;
using namespace mfem;

namespace {

struct CleanTrees {
	std::vector<std::filesystem::path> paths;
	~CleanTrees()
	{
		MPI_Barrier(MPI_COMM_WORLD);
		std::error_code ec;
		for (const auto& path : paths) {
			std::filesystem::remove_all(path, ec);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
};

void writeText(const std::filesystem::path& path, const std::string& text)
{
	std::filesystem::create_directories(path.parent_path());
	std::ofstream out(path, std::ios::trunc);
	ASSERT_TRUE(static_cast<bool>(out)) << path;
	out << text;
}

std::vector<std::string> readLines(const std::filesystem::path& path)
{
	std::ifstream in(path);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(in, line)) {
		lines.push_back(line);
	}
	return lines;
}

int countCheckpointDirs(const std::filesystem::path& root)
{
	int count = 0;
	if (!std::filesystem::exists(root)) {
		return 0;
	}
	for (const auto& entry : std::filesystem::directory_iterator(root)) {
		if (entry.is_directory()
			&& entry.path().filename().string().rfind("checkpoint_", 0) == 0) {
			++count;
		}
	}
	return count;
}

CheckpointManifest sampleManifest()
{
	CheckpointManifest meta;
	meta.time = 3.3;
	meta.dt = 0.1;
	meta.final_time = 10.0;
	meta.order = 2;
	meta.cycle = 12;
	meta.next_checkpoint_mark = 4;
	meta.elapsed_run_seconds = 240.0;
	ExporterCursor exporter;
	exporter.name = "view";
	exporter.save_count = 2;
	exporter.next_save_time = 4.0;
	exporter.dt_save = 2.0;
	exporter.initialized = true;
	exporter.finished = false;
	meta.exporters.push_back(exporter);
	MorCursor mor;
	mor.name = "mor";
	mor.save_count = 5;
	mor.next_save_time = 3.3;
	mor.dt_save = 0.5;
	mor.initialized = true;
	meta.mor.push_back(mor);
	return meta;
}

}  // namespace

class CheckpointTest : public ::testing::Test {
};

TEST_F(CheckpointTest, sha256OfEmptyAndAbc)
{
	const auto dir = std::filesystem::temp_directory_path()
		/ ("dgtd-ckpt-sha-" + std::to_string(Mpi::WorldRank()));
	CleanTrees clean{{dir}};

	const auto empty = dir / "empty";
	const auto abc = dir / "abc";
	writeText(empty, "");
	{
		std::ofstream out(abc, std::ios::binary | std::ios::trunc);
		out.write("abc", 3);
	}

	EXPECT_EQ(sha256File(empty),
		"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
	EXPECT_EQ(sha256File(abc),
		"ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}

TEST_F(CheckpointTest, roundTripStatePartitionAndSgbc)
{
	const std::string case_name = "_ckpt_test_roundtrip";
	const auto scratch = std::filesystem::temp_directory_path() / "dgtd-ckpt-roundtrip";
	CleanTrees clean{{
		std::filesystem::path(getSimulationCaseExportPath(case_name)),
		scratch
	}};

	if (Mpi::WorldRank() == 0) {
		writeText(scratch / "case.json", "{}\n");
		writeText(scratch / "mesh.msh", "mesh\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	const int rank = Mpi::WorldRank();
	const std::vector<double> state = {10.0 + rank, -2.5, 0.125 * (rank + 1)};
	SgbcRecord sgbc;
	sgbc.tag = 4 + rank;
	sgbc.node_a = 8;
	sgbc.node_b = -1;
	sgbc.values = {1.0, rank + 0.5};
	const std::vector<int> partition = {0, 1, 0, 2};

	writeCheckpointCollective(
		case_name,
		(scratch / "case.json").string(),
		(scratch / "mesh.msh").string(),
		sampleManifest(),
		partition.data(),
		static_cast<int>(partition.size()),
		state.data(),
		static_cast<int>(state.size()),
		{sgbc});

	const auto saved = findLatestCompleteCheckpoint(case_name);
	ASSERT_TRUE(saved.has_value());
	const auto manifest = readManifest(*saved / "manifest.json");
	EXPECT_EQ(manifest.world_size, Mpi::WorldSize());
	EXPECT_DOUBLE_EQ(manifest.time, 3.3);
	EXPECT_DOUBLE_EQ(manifest.elapsed_run_seconds, 240.0);
	EXPECT_EQ(manifest.cycle, 12);
	EXPECT_EQ(manifest.next_checkpoint_mark, 4);
	ASSERT_EQ(static_cast<int>(manifest.state_sizes.size()), Mpi::WorldSize());
	EXPECT_EQ(manifest.state_sizes[static_cast<std::size_t>(rank)], 3);
	EXPECT_EQ(manifest.sgbc_counts[static_cast<std::size_t>(rank)], 1);
	ASSERT_EQ(static_cast<int>(manifest.state_sha256.size()), Mpi::WorldSize());
	EXPECT_EQ(manifest.state_sha256[static_cast<std::size_t>(rank)].size(), 64u);
	ASSERT_EQ(manifest.exporters.size(), 1u);
	EXPECT_EQ(manifest.exporters[0].name, "view");
	EXPECT_FALSE(manifest.exporters[0].finished);
	ASSERT_EQ(manifest.mor.size(), 1u);
	EXPECT_EQ(manifest.mor[0].save_count, 5);
	EXPECT_EQ(readPartition(*saved / "partition.bin"), partition);
	EXPECT_EQ(sha256File(*saved / "input.json"), manifest.json_sha256);

	const auto state_path = *saved / ("state.rank" + std::to_string(rank) + ".bin");
	EXPECT_EQ(sha256File(state_path), manifest.state_sha256[static_cast<std::size_t>(rank)]);

	std::vector<double> loaded(state.size());
	std::vector<SgbcRecord> loaded_sgbc;
	readRankState(state_path, loaded.data(), static_cast<int>(loaded.size()), loaded_sgbc);
	EXPECT_EQ(loaded, state);
	ASSERT_EQ(loaded_sgbc.size(), 1u);
	EXPECT_EQ(loaded_sgbc[0].tag, sgbc.tag);
	EXPECT_EQ(loaded_sgbc[0].node_a, 8);
	EXPECT_EQ(loaded_sgbc[0].node_b, -1);
	EXPECT_EQ(loaded_sgbc[0].values, sgbc.values);

	const auto corrupt = scratch / ("rank" + std::to_string(rank) + "-bad-magic.bin");
	std::filesystem::copy_file(state_path, corrupt);
	{
		std::fstream flipped(corrupt, std::ios::binary | std::ios::in | std::ios::out);
		ASSERT_TRUE(static_cast<bool>(flipped));
		char byte = 0;
		flipped.read(&byte, 1);
		flipped.seekp(0);
		byte = static_cast<char>(byte ^ 0x1);
		flipped.write(&byte, 1);
	}
	EXPECT_THROW(
		readRankState(corrupt, loaded.data(), static_cast<int>(loaded.size()), loaded_sgbc),
		std::runtime_error);

	const auto tampered = scratch / ("rank" + std::to_string(rank) + "-tampered.bin");
	std::filesystem::copy_file(state_path, tampered);
	{
		std::fstream flipped(tampered, std::ios::binary | std::ios::in | std::ios::out);
		ASSERT_TRUE(static_cast<bool>(flipped));
		flipped.seekp(20);
		char byte = 0;
		flipped.read(&byte, 1);
		flipped.seekp(20);
		byte = static_cast<char>(byte ^ 0x1);
		flipped.write(&byte, 1);
	}
	EXPECT_NE(sha256File(tampered), manifest.state_sha256[static_cast<std::size_t>(rank)]);
}

TEST_F(CheckpointTest, completeMarkerKeepsNewestAndDropsTheRest)
{
	const std::string case_name = "_ckpt_test_generations";
	const auto root = checkpointRoot(case_name);
	const auto scratch = std::filesystem::temp_directory_path() / "dgtd-ckpt-generations";
	CleanTrees clean{{std::filesystem::path(getSimulationCaseExportPath(case_name)), scratch}};

	if (Mpi::WorldRank() == 0) {
		writeText(scratch / "case.json", "{}\n");
		writeText(scratch / "mesh.msh", "mesh\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	const double state = 1.0;
	const int partition = 0;
	auto write_one = [&]() {
		writeCheckpointCollective(
			case_name,
			(scratch / "case.json").string(),
			(scratch / "mesh.msh").string(),
			CheckpointManifest{},
			&partition,
			1,
			&state,
			1,
			{});
	};

	write_one();
	const auto first = findLatestCompleteCheckpoint(case_name);
	ASSERT_TRUE(first.has_value());

	if (Mpi::WorldRank() == 0) {
		const auto abandoned = root / "checkpoint_000050";
		std::filesystem::create_directories(abandoned);
		writeText(abandoned / "manifest.json", "{}\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	EXPECT_EQ(findLatestCompleteCheckpoint(case_name), first);

	write_one();
	const auto newest = findLatestCompleteCheckpoint(case_name);
	ASSERT_TRUE(newest.has_value());
	EXPECT_NE(*newest, *first);
	EXPECT_EQ(newest->filename(), "checkpoint_000051");
	EXPECT_TRUE(std::filesystem::exists(*newest / "COMPLETE"));
	EXPECT_FALSE(std::filesystem::exists(root / "checkpoint_000050"));
	EXPECT_FALSE(std::filesystem::exists(*first));
	EXPECT_EQ(countCheckpointDirs(root), 1);
}

TEST_F(CheckpointTest, manifestWithoutChecksumStillLoads)
{
	const auto path = std::filesystem::temp_directory_path()
		/ ("dgtd-ckpt-manifest-" + std::to_string(Mpi::WorldRank()) + ".json");
	CleanTrees clean{{path}};
	writeText(path, R"({
  "format_version": 1,
  "json_sha256": "abc",
  "mesh_sha256": "def",
  "world_size": 1,
  "time": 1.5,
  "dt": 0.1,
  "final_time": 4.0,
  "order": 1,
  "cycle": 3,
  "next_checkpoint_mark": 2,
  "state_sizes": [6],
  "sgbc_counts": [0],
  "exporters": [],
  "mor": []
}
)");

	const auto manifest = readManifest(path);
	EXPECT_TRUE(manifest.state_sha256.empty());
	EXPECT_DOUBLE_EQ(manifest.time, 1.5);
	EXPECT_DOUBLE_EQ(manifest.elapsed_run_seconds, 0.0);
	EXPECT_EQ(manifest.state_sizes, std::vector<int>({6}));
}

TEST_F(CheckpointTest, foreignRunModeIsReported)
{
	const std::string case_name = "_ckpt_test_foreign";
	const auto foreign_root = std::filesystem::path("exports/SimulationData");
	const auto mpi99 = foreign_root / "mpi-99" / case_name;
	const auto cuda1 = foreign_root / "cuda-1" / case_name;
	CleanTrees clean{{mpi99, cuda1}};
	if (Mpi::WorldRank() == 0) {
		std::error_code ec;
		std::filesystem::remove_all(mpi99, ec);
		std::filesystem::remove_all(cuda1, ec);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	const std::string foreign_mode = (getRunModeTag() == "mpi-99") ? "cuda-1" : "mpi-99";
	const auto foreign = foreign_root / foreign_mode / case_name;

	const auto checkpoint = foreign / "Checkpoints" / "checkpoint_000003";
	if (Mpi::WorldRank() == 0) {
		std::filesystem::create_directories(checkpoint);
		writeText(checkpoint / "COMPLETE", "ok\n");
		writeText(checkpoint / "manifest.json", R"({
  "format_version": 1,
  "json_sha256": "abc",
  "mesh_sha256": "def",
  "world_size": 6,
  "time": 1.0,
  "dt": 0.1,
  "final_time": 10.0,
  "order": 2,
  "cycle": 1,
  "next_checkpoint_mark": 1,
  "state_sizes": [1, 1, 1, 1, 1, 1],
  "sgbc_counts": [0, 0, 0, 0, 0, 0],
  "exporters": [],
  "mor": []
}
)");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	const auto found = findForeignCompleteCheckpoint(case_name);
	ASSERT_TRUE(found.has_value());
	EXPECT_EQ(*found, checkpoint);
	EXPECT_FALSE(findLatestCompleteCheckpoint(case_name).has_value());

	const auto message = foreignCheckpointMessage(*found);
	EXPECT_NE(message.find(checkpoint.string()), std::string::npos);
	EXPECT_NE(message.find("(6 ranks)"), std::string::npos);
	EXPECT_NE(message.find("This job is " + getRunModeTag()), std::string::npos);
}

TEST_F(CheckpointTest, percentMarksSnapAtFinalTime)
{
	const auto at_33 = checkpointSchedule(33.0, 0.1, 100.0, 10.0, 3);
	EXPECT_TRUE(at_33.due);
	EXPECT_EQ(at_33.next_mark, 4);

	const auto already_saved = checkpointSchedule(33.0, 0.1, 100.0, 10.0, 4);
	EXPECT_FALSE(already_saved.due);
	EXPECT_EQ(already_saved.next_mark, 4);

	const auto missed = checkpointSchedule(33.0, 0.1, 100.0, 10.0, 1);
	EXPECT_TRUE(missed.due);
	EXPECT_EQ(missed.next_mark, 4);

	const auto disabled = checkpointSchedule(33.0, 0.1, 100.0, 0.0, 1);
	EXPECT_FALSE(disabled.due);
	EXPECT_EQ(disabled.next_mark, 1);

	const auto at_final = checkpointSchedule(1.0, 0.1, 1.0, 10.0, 10);
	EXPECT_TRUE(at_final.due);
	EXPECT_EQ(at_final.next_mark, 11);

	const auto inside_snap = checkpointSchedule(1.0 - 5e-9, 1.0, 1.0, 20.0, 5);
	EXPECT_TRUE(inside_snap.due);
	EXPECT_EQ(inside_snap.next_mark, 6);

	const auto short_of_final = checkpointSchedule(1.0 - 1e-6, 0.1, 1.0, 10.0, 10);
	EXPECT_FALSE(short_of_final.due);
	EXPECT_EQ(short_of_final.next_mark, 10);
}

TEST_F(CheckpointTest, probeSampleCountAndPvdCycle)
{
	EXPECT_EQ(expectedProbeSamples(0, 10, false), 0);
	EXPECT_EQ(expectedProbeSamples(1, 10, false), 1);
	EXPECT_EQ(expectedProbeSamples(11, 10, false), 2);
	EXPECT_EQ(expectedProbeSamples(11, 10, true), 2);
	EXPECT_EQ(expectedProbeSamples(12, 10, true), 3);
	EXPECT_EQ(expectedProbeSamples(12, 0, false), 0);
	EXPECT_EQ(expectedProbeSamples(12, 0, true), 1);

	EXPECT_EQ(cycleInPvdLine("Cycle000012/data.pvtu"), 12);
	EXPECT_EQ(cycleInPvdLine("<DataSet file=\"Cycle000012/data.pvtu\"/>"), 12);
	EXPECT_EQ(cycleInPvdLine("Cycle12abc"), 12);
	EXPECT_EQ(cycleInPvdLine("<?xml version=\"1.0\"?>"), -1);
	EXPECT_EQ(cycleInPvdLine("Cycle"), -1);
}

TEST_F(CheckpointTest, jsonMayChangeOnlyCheckpointPercent)
{
	const auto dir = std::filesystem::temp_directory_path()
		/ ("dgtd-ckpt-json-" + std::to_string(Mpi::WorldRank()));
	CleanTrees clean{{dir}};

	const auto saved = dir / "saved.json";
	const auto percent_only = dir / "percent.json";
	const auto other = dir / "other.json";
	writeText(saved, R"({"solver_options":{"final_time":1.0,"checkpoint_percent":20},"model":{"filename":"a.msh"}})");
	writeText(percent_only, R"({"solver_options":{"checkpoint_percent":50,"final_time":1.0},"model":{"filename":"a.msh"}})");
	writeText(other, R"({"solver_options":{"final_time":2.0,"checkpoint_percent":20},"model":{"filename":"a.msh"}})");

	EXPECT_TRUE(sameCaseExceptCheckpointPercent(percent_only, saved));
	EXPECT_FALSE(sameCaseExceptCheckpointPercent(other, saved));
	EXPECT_FALSE(sameCaseExceptCheckpointPercent(dir / "missing.json", saved));
}

TEST_F(CheckpointTest, trimDropsProbeTailAndNewerParaviewCycles)
{
	if (Mpi::WorldSize() != 1) {
		return;
	}

	const std::string case_name = "_ckpt_test_trim";
	const std::string probe_name = "_ckpt_test_pvd";
	CleanTrees clean{{
		std::filesystem::path(getSimulationCaseExportPath(case_name)),
		std::filesystem::path("exports/ParaView") / getRunModeTag() / probe_name
	}};

	Mesh smesh = Mesh::MakeCartesian1D(5, 1.0);
	ParMesh mesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec(2, 1, BasisType::GaussLobatto);
	ParFiniteElementSpace fes(&mesh, &fec);
	Fields<ParFiniteElementSpace, ParGridFunction> fields(fes);

	PointProbe point({0.5}, 10);
	point.setProbeID(0);
	ExporterProbe exporter;
	exporter.name = probe_name;
	exporter.save_every = 1.0;

	Probes probes;
	probes.pointProbes = {point};
	probes.exporterProbes = {exporter};

	SolverOptions opts;
	opts.final_time = 1.0;
	ProbesManager manager(probes, fes, fields, opts);
	manager.setCaseName(case_name);

	const auto dat = std::filesystem::path(getSimulationCaseExportPath(case_name))
		/ "PointProbes" / "PointProbe0.dat";
	writeText(dat,
		"h0\nh1\nh2\nh3\n"
		"d0\nd1\nd2\nd3\nd4\nd5\n");

	const auto paraview = std::filesystem::path("exports/ParaView") / getRunModeTag() / probe_name;
	std::filesystem::create_directories(paraview / "Cycle000000");
	std::filesystem::create_directories(paraview / "Cycle000001");
	std::filesystem::create_directories(paraview / "Cycle000002");
	writeText(paraview / (probe_name + ".pvd"),
		"<?xml version=\"1.0\"?>\n"
		"<DataSet file=\"Cycle000000/data.pvtu\"/>\n"
		"<DataSet file=\"Cycle000001/data.pvtu\"/>\n"
		"<DataSet file=\"Cycle000002/data.pvtu\"/>\n");

	ExporterCursor cursor;
	cursor.name = probe_name;
	cursor.save_count = 2;
	cursor.next_save_time = 2.0;
	cursor.dt_save = 1.0;
	cursor.initialized = true;
	manager.restoreCheckpointCursors(12, {cursor}, {});
	manager.trimProbeOutput(1.0);

	const auto kept = readLines(dat);
	ASSERT_EQ(kept.size(), 7u);
	EXPECT_EQ(kept.back(), "d2");

	EXPECT_TRUE(std::filesystem::exists(paraview / "Cycle000000"));
	EXPECT_TRUE(std::filesystem::exists(paraview / "Cycle000001"));
	EXPECT_FALSE(std::filesystem::exists(paraview / "Cycle000002"));
	const auto pvd = readLines(paraview / (probe_name + ".pvd"));
	ASSERT_EQ(pvd.size(), 3u);
	EXPECT_EQ(cycleInPvdLine(pvd[1]), 0);
	EXPECT_EQ(cycleInPvdLine(pvd[2]), 1);
}
