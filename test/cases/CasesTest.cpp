#include "driver/driver.h"

#include "TestUtils.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;

using namespace mfem;
using namespace maxwell;
using namespace maxwell::driver;

class CasesTest : public ::testing::Test {

};

TEST_F(CasesTest, 1D_PEC_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("1D_PEC"));
    case_data["solver_options"]["upwind_alpha"] = 0.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_PEC"), true) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}

}

TEST_F(CasesTest, 1D_PEC_Upwind)
{

	std::string case_name{ "1D_PEC" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}

}


TEST_F(CasesTest, 1D_SMA_Upwind)
{
	auto case_data = parseJSONfile(maxwellCase("1D_PEC"));
    case_data["model"]["boundaries"][0]["type"] = "SMA";
	case_data["solver_options"]["final_time"] = 1.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_PEC"), true) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 1e-2);

}


TEST_F(CasesTest, 1D_MultiCondition_Upwind)
{
	SGBCProperties sbcp;
	SGBCLayer layer(Material(20.0, 1.0, 20.0), 2.0);
	layer.order = 4;
	layer.num_of_segments = 10;
	sbcp.layers.push_back(layer);
	SGBCBoundaries sbcp_bdrs;
	sbcp_bdrs.second.isOn = true;
	sbcp_bdrs.second.bdrCond = BdrCond::PEC;
	sbcp.sgbc_bdr_info = sbcp_bdrs;

	auto total_segments = sbcp.totalSegments();
	auto total_width = sbcp.totalWidth();
	auto mesh = mfem::Mesh::MakeCartesian1D(total_segments + 2, total_width + 2 * total_width / total_segments);
    mesh.AddBdrPoint(1, 3);
    mesh.AddBdrPoint(mesh.GetNV() - 2, 4);
    mesh.SetAttribute(0, 2);
    mesh.SetAttribute(mesh.GetNE() - 1, 2);
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); 
    mesh.FinalizeMesh();
    int* partitioning = mesh.GeneratePartitioning(1);

    Material vacuum = buildVacuumMaterial();
    GeomTagToMaterial geom_tag_sgbc_mat{{1, vacuum}, {2, vacuum}};
    GeomTagToInteriorBoundary gt2ib;
    if (sbcp_bdrs.first.isOn){
        gt2ib[3] = sbcp_bdrs.first.bdrCond;
    }
    if (sbcp_bdrs.second.isOn){
        gt2ib[4] = sbcp_bdrs.second.bdrCond;
    }
	GeomTagToBoundary gt2b;
    gt2b[1] = BdrCond::SMA;
    gt2b[2] = BdrCond::SMA;
    GeomTagToBoundaryInfo gtbdr(gt2b, gt2ib);
    Model model = Model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), gtbdr, partitioning);

	Probes probes;
    probes.exporterProbes.resize(1);
    ExporterProbe ep;
    ep.name = "1D_MultiCondition_Upwind";
    ep.visSteps = 100;
    probes.exporterProbes.at(0) = ep;
    Sources sources;
	mfem::Vector gaussianCenter(1);
	gaussianCenter = 0.5;
	double spread = 0.1;
	Gaussian gauss{ spread, gaussianCenter, 1 };
	sources.add(std::make_unique<InitialField>(gauss, E, mfem::Vector({0.0,1.0,0.0}), gaussianCenter));
	sources.add(std::make_unique<InitialField>(gauss, H, mfem::Vector({0.0,0.0,1.0}), gaussianCenter));
    SolverOptions opts;
    opts.setOrder(sbcp.maxOrder());
    opts.setUpwindAlpha(1.0);
    opts.setODEType(ode_type::RK4); 
	opts.setFinalTime(10.0);
	opts.setTimeStep(1e-3);
    opts.setExportEO(true);
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    auto solver = maxwell::Solver(model, probes, sources, opts);

	solver.run();

}


TEST_F(CasesTest, 1D_TFSF_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("1D_TFSF"));
    case_data["solver_options"]["upwind_alpha"] = 0.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_TFSF"), true) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 4.2 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 10.2 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 7.2 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}

TEST_F(CasesTest, 1D_TFSF_Upwind)
{
	std::string case_name{ "1D_TFSF" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 4.2 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 10.2 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 7.2 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}

TEST_F(CasesTest, 1D_TFSF_Upwind_Negative_K)
{
	auto case_data = parseJSONfile(maxwellCase("1D_TFSF"));
	case_data["probes"]["point"][0]["position"] = {2.99999};
	case_data["probes"]["point"][2]["position"] = { -2.99999 };
	case_data["sources"][0]["propagation"] = { -1.0, 0.0, 0.0 };
	case_data["sources"][0]["magnitude"]["mean"] = { 5.0, 0.0, 0.0 };
	auto solver{ buildSolver(case_data, maxwellCase("1D_TFSF"), true) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}

namespace {

struct SkrfReferenceRow {
	double freq_hz{0.0};
	double s11_db{0.0};
	double s21_db{0.0};
	bool has_s21{false};
};

std::vector<SkrfReferenceRow> loadSkrfReferenceCsv(const std::string& csv_path)
{
	std::ifstream in(csv_path);
	if (!in.is_open()) {
		throw std::runtime_error("Cannot open scikit-rf reference CSV: " + csv_path);
	}

	std::vector<SkrfReferenceRow> rows;
	std::string line;
	// Skip header
	if (!std::getline(in, line)) {
		throw std::runtime_error("Empty scikit-rf reference CSV: " + csv_path);
	}

	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		std::stringstream ss(line);
		std::string token;
		SkrfReferenceRow row;
		if (!std::getline(ss, token, ',')) {
			continue;
		}
		row.freq_hz = std::stod(token);
		if (!std::getline(ss, token, ',')) {
			throw std::runtime_error("Missing s11_db in CSV: " + csv_path);
		}
		row.s11_db = std::stod(token);
		if (std::getline(ss, token, ',')) {
			row.s21_db = std::stod(token);
			row.has_s21 = true;
		}
		rows.push_back(row);
	}
	return rows;
}

std::vector<std::complex<double>> computeProbeDftEy(
	const FieldMovies& movie,
	const std::vector<SkrfReferenceRow>& ref_rows,
	double max_code_time = std::numeric_limits<double>::infinity())
{
	constexpr double c_si = physicalConstants::speedOfLight_SI;
	std::vector<std::complex<double>> spectrum(ref_rows.size(), {0.0, 0.0});

	for (const auto& [t_code, f] : movie) {
		if (t_code >= max_code_time || f.Ey == 0.0) {
			continue;
		}
		const double t_si = t_code / c_si;
		for (std::size_t i = 0; i < ref_rows.size(); ++i) {
			const double angle = -2.0 * M_PI * ref_rows[i].freq_hz * t_si;
			spectrum[i] += f.Ey * std::complex<double>(std::cos(angle), std::sin(angle));
		}
	}

	if (Mpi::WorldSize() > 1) {
		struct RankCount {
			int count;
			int rank;
		} local{static_cast<int>(movie.size()), Mpi::WorldRank()}, chosen{0, 0};
		MPI_Allreduce(&local, &chosen, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
		if (chosen.count > 0 && !spectrum.empty()) {
			MPI_Bcast(
				reinterpret_cast<double*>(spectrum.data()),
				static_cast<int>(2 * spectrum.size()),
				MPI_DOUBLE,
				chosen.rank,
				MPI_COMM_WORLD);
		}
	}

	return spectrum;
}

double computeProbeMaxAbsEy(const FieldMovies& movie)
{
	double local_max = 0.0;
	for (const auto& [t_code, f] : movie) {
		(void)t_code;
		local_max = std::max(local_max, std::abs(f.Ey));
	}
	double global_max = local_max;
	if (Mpi::WorldSize() > 1) {
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	return global_max;
}

} // namespace

TEST_F(CasesTest, 1D_SGBC_Multilayer_PEC_Backed_RAM_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Multilayer"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Multilayer_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	// Probe 0: reflected region (x = 0.10 m)
	// Probe 1: incident region (x = 0.50 m), time-gated to t_code < 3.5 m before reflection arrives
	// Probe 2: behind PEC-backed SGBC (x = 3.45 m)
	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	EXPECT_NEAR(0.0, computeProbeMaxAbsEy(trn_movie), 1e-12);

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 3.5);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);

	double sum_abs_err_db = 0.0;
	double max_abs_err_db = 0.0;
	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double err_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		sum_abs_err_db += err_db;
		max_abs_err_db = std::max(max_abs_err_db, err_db);
		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.10);
	}

	const double mean_abs_err_db = sum_abs_err_db / static_cast<double>(ref_rows.size());
	// Paper Section IV-A: mean error ~0.035 dB, max error ~0.075 dB over 4.5-7.4 GHz
	EXPECT_LT(mean_abs_err_db, 0.05);
	EXPECT_LT(max_abs_err_db, 0.08);
}

TEST_F(CasesTest, 1D_SGBC_Resonant_Passthrough_Slab_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Resonant"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Resonant_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	// Probe 0: reflected region (x = 0.10 m)
	// Probe 1: incident region (x = 0.50 m), time-gated to t_code < 3.5 m before reflection arrives
	// Probe 2: transmitted region behind SGBC (x = 3.45 m)
	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 3.5);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);
	const auto F_trn = computeProbeDftEy(trn_movie, ref_rows);

	double sum_s11_err_db = 0.0;
	double max_s11_err_db = 0.0;
	double sum_s21_err_db = 0.0;
	double max_s21_err_db = 0.0;

	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_TRUE(ref_rows[i].has_s21);
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double sim_s21_db = 20.0 * std::log10(std::abs(F_trn[i]) / std::abs(F_inc[i]));

		const double err_s11_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		const double err_s21_db = std::abs(sim_s21_db - ref_rows[i].s21_db);

		sum_s11_err_db += err_s11_db;
		max_s11_err_db = std::max(max_s11_err_db, err_s11_db);
		sum_s21_err_db += err_s21_db;
		max_s21_err_db = std::max(max_s21_err_db, err_s21_db);

		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.03);
		EXPECT_NEAR(ref_rows[i].s21_db, sim_s21_db, 0.01);
	}

	const double mean_s11_err_db = sum_s11_err_db / static_cast<double>(ref_rows.size());
	const double mean_s21_err_db = sum_s21_err_db / static_cast<double>(ref_rows.size());
	// Paper Section IV-B: S11 mean error ~0.011 dB, max ~0.015 dB; S21 within noise floor over 5.0-7.0 GHz
	EXPECT_LT(mean_s11_err_db, 0.02);
	EXPECT_LT(max_s11_err_db, 0.025);
	EXPECT_LT(mean_s21_err_db, 0.005);
	EXPECT_LT(max_s21_err_db, 0.01);
}

TEST_F(CasesTest, 1D_SGBC_Through_Passthrough_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Through"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Through_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	// Probe 0: reflected (x = 0.15 m), Probe 1: incident (x = 0.35 m, gated t_code < 1.5 m),
	// Probe 2: transmitted (x = 1.45 m)
	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 1.5);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);
	const auto F_trn = computeProbeDftEy(trn_movie, ref_rows);

	double sum_s11_err_db = 0.0;
	double max_s11_err_db = 0.0;
	double sum_s21_err_db = 0.0;
	double max_s21_err_db = 0.0;

	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_TRUE(ref_rows[i].has_s21);
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double sim_s21_db = 20.0 * std::log10(std::abs(F_trn[i]) / std::abs(F_inc[i]));

		const double err_s11_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		const double err_s21_db = std::abs(sim_s21_db - ref_rows[i].s21_db);

		sum_s11_err_db += err_s11_db;
		max_s11_err_db = std::max(max_s11_err_db, err_s11_db);
		sum_s21_err_db += err_s21_db;
		max_s21_err_db = std::max(max_s21_err_db, err_s21_db);

		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.06);
		EXPECT_NEAR(ref_rows[i].s21_db, sim_s21_db, 0.03);
	}

	EXPECT_LT(sum_s11_err_db / static_cast<double>(ref_rows.size()), 0.03);
	EXPECT_LT(max_s11_err_db, 0.06);
	EXPECT_LT(sum_s21_err_db / static_cast<double>(ref_rows.size()), 0.01);
	EXPECT_LT(max_s21_err_db, 0.03);
}

TEST_F(CasesTest, 1D_SGBC_Through_PEC_Backed_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Through"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}
	case_data["model"]["boundaries"][1]["material"]["relative_permittivity"] = 4.0;
	case_data["model"]["boundaries"][1]["material"]["bulk_conductivity"] = 10.0;
	case_data["model"]["boundaries"][1]["sgbc_boundaries"]["right"] = "PEC";

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Through_PEC_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	EXPECT_NEAR(0.0, computeProbeMaxAbsEy(trn_movie), 1e-12);

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 1.5);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);

	double sum_s11_err_db = 0.0;
	double max_s11_err_db = 0.0;
	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double err_s11_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		sum_s11_err_db += err_s11_db;
		max_s11_err_db = std::max(max_s11_err_db, err_s11_db);
		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.08);
	}

	EXPECT_LT(sum_s11_err_db / static_cast<double>(ref_rows.size()), 0.05);
	EXPECT_LT(max_s11_err_db, 0.08);
}

TEST_F(CasesTest, 1D_SGBC_Multiple_Passthrough_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Multiple"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Multiple_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	// Probe 0: reflected (x = 0.10 m)
	// Probe 1: incident (x = 0.35 m, gated t_code < 1.45 m before first SGBC reflection at x = 1.00 m arrives)
	// Probe 2: transmitted after all 3 SGBCs (x = 2.25 m)
	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 1.45);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);
	const auto F_trn = computeProbeDftEy(trn_movie, ref_rows);

	double sum_s11_err_db = 0.0;
	double max_s11_err_db = 0.0;
	double sum_s21_err_db = 0.0;
	double max_s21_err_db = 0.0;

	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_TRUE(ref_rows[i].has_s21);
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double sim_s21_db = 20.0 * std::log10(std::abs(F_trn[i]) / std::abs(F_inc[i]));

		const double err_s11_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		const double err_s21_db = std::abs(sim_s21_db - ref_rows[i].s21_db);

		sum_s11_err_db += err_s11_db;
		max_s11_err_db = std::max(max_s11_err_db, err_s11_db);
		sum_s21_err_db += err_s21_db;
		max_s21_err_db = std::max(max_s21_err_db, err_s21_db);

		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.12);
		EXPECT_NEAR(ref_rows[i].s21_db, sim_s21_db, 0.10);
	}

	EXPECT_LT(sum_s11_err_db / static_cast<double>(ref_rows.size()), 0.05);
	EXPECT_LT(max_s11_err_db, 0.12);
	EXPECT_LT(sum_s21_err_db / static_cast<double>(ref_rows.size()), 0.04);
	EXPECT_LT(max_s21_err_db, 0.10);
}

TEST_F(CasesTest, 1D_SGBC_Multiple_PEC_Backed_vs_Skrf)
{
	const std::string case_name{"1D_SGBC_Multiple"};
	auto case_data = parseJSONfile(maxwellCase(case_name));
	case_data["probes"].erase("exporter");
	for (auto& probe_json : case_data["probes"]["point"]) {
		probe_json["write"] = false;
	}
	// Terminate the third SGBC boundary (tag 5 at x = 2.00 m) with a PEC backing
	case_data["model"]["boundaries"][3]["sgbc_boundaries"]["right"] = "PEC";

	auto solver{ buildSolver(case_data, maxwellCase(case_name), true) };
	solver.run();

	const auto ref_rows = loadSkrfReferenceCsv(
		testDataFolder() + "sgbcReferences/1D_SGBC_Multiple_PEC_skrf.csv");
	ASSERT_FALSE(ref_rows.empty());

	const auto& ref_movie = solver.getPointProbe(0).getFieldMovies();
	const auto& inc_movie = solver.getPointProbe(1).getFieldMovies();
	const auto& trn_movie = solver.getPointProbe(2).getFieldMovies();

	EXPECT_NEAR(0.0, computeProbeMaxAbsEy(trn_movie), 1e-12);

	const auto F_inc = computeProbeDftEy(inc_movie, ref_rows, 1.45);
	const auto F_ref = computeProbeDftEy(ref_movie, ref_rows);

	double sum_s11_err_db = 0.0;
	double max_s11_err_db = 0.0;
	for (std::size_t i = 0; i < ref_rows.size(); ++i) {
		ASSERT_GT(std::abs(F_inc[i]), 1e-12);
		const double sim_s11_db = 20.0 * std::log10(std::abs(F_ref[i]) / std::abs(F_inc[i]));
		const double err_s11_db = std::abs(sim_s11_db - ref_rows[i].s11_db);
		sum_s11_err_db += err_s11_db;
		max_s11_err_db = std::max(max_s11_err_db, err_s11_db);
		EXPECT_NEAR(ref_rows[i].s11_db, sim_s11_db, 0.18);
	}

	EXPECT_LT(sum_s11_err_db / static_cast<double>(ref_rows.size()), 0.05);
	EXPECT_LT(max_s11_err_db, 0.18);
}