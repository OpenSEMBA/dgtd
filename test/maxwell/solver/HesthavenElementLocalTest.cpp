#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

#include "driver/driver.h"

#include "TestUtils.h"
#include "components/DebyeProperties.h"
#include "components/LorentzProperties.h"
#include "components/Model.h"
#include "components/PMLProperties.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolution.h"
#include "math/Function.h"
#include "solver/SourcesManager.h"

using namespace mfem;
using namespace maxwell;

namespace {

GeomTagToBoundary pecBoundary(const Mesh& mesh)
{
	GeomTagToBoundary bdr;
	for (int i = 0; i < mesh.GetNBE(); ++i) {
		bdr[mesh.GetBdrAttribute(i)] = BdrCond::PEC;
	}
	return bdr;
}

GeomTagToBoundary smaBoundary(const Mesh& mesh)
{
	GeomTagToBoundary bdr;
	for (int i = 0; i < mesh.GetNBE(); ++i) {
		bdr[mesh.GetBdrAttribute(i)] = BdrCond::SMA;
	}
	return bdr;
}

void elementCentroid(const Mesh& mesh, int elem, double& x, double& y)
{
	Array<int> verts;
	mesh.GetElementVertices(elem, verts);
	x = 0.0;
	y = 0.0;
	for (int v = 0; v < verts.Size(); ++v) {
		const double* coord = mesh.GetVertex(verts[v]);
		x += coord[0];
		y += coord[1];
	}
	x /= static_cast<double>(verts.Size());
	y /= static_cast<double>(verts.Size());
}

void tagRightHalf(Mesh& mesh, int right_attribute)
{
	for (int elem = 0; elem < mesh.GetNE(); ++elem) {
		double x = 0.0;
		double y = 0.0;
		elementCentroid(mesh, elem, x, y);
		if (x >= 0.5) {
			mesh.SetAttribute(elem, right_attribute);
		}
	}
	mesh.SetAttributes();
}

struct PairRun {
	std::unique_ptr<Model> model;
	std::unique_ptr<DG_FECollection> fec;
	std::unique_ptr<ParFiniteElementSpace> fes;
	Sources sources;
	std::unique_ptr<Fields<ParFiniteElementSpace, ParGridFunction>> fields;
	std::unique_ptr<SourcesManager> sources_manager;
	EvolutionOptions opts;
	Probes probes;
	std::unique_ptr<GlobalEvolution> global;
	std::unique_ptr<HesthavenEvolution> hesthaven;

	void open(int order, double alpha = 1.0)
	{
		const int dim = model->getMesh().Dimension();
		fec = std::make_unique<DG_FECollection>(order, dim, BasisType::GaussLobatto);
		fes = std::make_unique<ParFiniteElementSpace>(&model->getMesh(), fec.get());
		const int n_aux = computePMLAuxSize(model->getPMLProperties(), fes->GetNDofs(), dim)
			+ model->debyeAuxSize(fes->GetNDofs())
			+ model->lorentzAuxSize(fes->GetNDofs());
		fields = std::make_unique<Fields<ParFiniteElementSpace, ParGridFunction>>(*fes, n_aux);
		sources_manager = std::make_unique<SourcesManager>(sources, *fes, *fields);
		opts.order = order;
		opts.alpha = alpha;
		global = std::make_unique<GlobalEvolution>(
			*fes, *model, *sources_manager, opts, probes, 1.0);
		hesthaven = std::make_unique<HesthavenEvolution>(
			*fes, *model, *sources_manager, opts);
	}
};

double maxAbsDiff(const Vector& a, const Vector& b)
{
	const double* ah = a.HostRead();
	const double* bh = b.HostRead();
	double diff = 0.0;
	for (int i = 0; i < a.Size(); ++i) {
		diff = std::max(diff, std::abs(ah[i] - bh[i]));
	}
	return diff;
}

void expectFullStateAgrees(PairRun& run, const char* name, bool aux_nonzero)
{
	const int height = run.global->Height();
	ASSERT_EQ(run.hesthaven->Height(), height) << name;
	Vector x(height);
	Vector y_global(height);
	Vector y_hest(height);
	for (int i = 0; i < height; ++i) {
		x[i] = 0.1 * std::sin(0.17 * static_cast<double>(i + 1));
	}
	auto expectAux = [&](const Vector& y) {
		if (!aux_nonzero) {
			return;
		}
		const int eh = 6 * run.fes->GetNDofs();
		const double* data = y.HostRead();
		double aux = 0.0;
		for (int i = eh; i < height; ++i) {
			aux += data[i] * data[i];
		}
		EXPECT_GT(aux, 1e-16) << name << " auxiliary residual is zero";
	};
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (mfem::Device::Allows(mfem::Backend::CUDA)) {
		// Host reads of GlobalEvolution's CUDA residual drop the auxiliary
		// blocks. The CPU Hesthaven residual matches Global on the host, so the
		// device kernel is checked against that CPU residual.
		Vector y_cpu(height);
		run.hesthaven->benchmarkMultCPU(x, y_cpu);
		run.hesthaven->Mult(x, y_hest);
		const double gpu_diff = maxAbsDiff(y_cpu, y_hest);
		EXPECT_LT(gpu_diff, 1e-8) << name << " max |y_cpu - y_gpu| = " << gpu_diff;
		expectAux(y_cpu);
		return;
	}
#endif
	run.global->Mult(x, y_global);
	run.hesthaven->Mult(x, y_hest);
	const double diff = maxAbsDiff(y_global, y_hest);
	EXPECT_LT(diff, 1e-8) << name << " max |y_global - y_hest| = " << diff;
	expectAux(y_global);
}

} // namespace

void expectImplicitAgrees(PairRun& run, const char* name);

TEST(HesthavenElementLocal, conductive_dielectric_matches_global)
{
	Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
	tagRightHalf(serial, 2);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	mats.emplace(2, Material(2.0, 1.5, 0.2));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	run.open(2);
	expectFullStateAgrees(run, "conductive dielectric", false);
}

TEST(HesthavenElementLocal, debye_matches_global)
{
	const double eps_inf = 2.0;
	Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = 6.0;
	pole.tau_solver = 1.0;
	run.model->setDebyeProperties({pole});
	run.open(2);
	expectFullStateAgrees(run, "debye", true);
}

TEST(HesthavenElementLocal, lorentz_matches_global)
{
	const double eps_inf = 2.0;
	Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.omega_p = 2.0;
	pole.omega_1 = 1.0;
	pole.gamma = 0.5;
	run.model->setLorentzProperties({pole});
	run.open(2);
	expectFullStateAgrees(run, "lorentz", true);
	expectImplicitAgrees(run, "implicit lorentz");
}

TEST(HesthavenElementLocal, profiled_drude_matches_global)
{
	Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = 1.0;
	pole.omega_1 = 0.0;
	pole.gamma = 0.4;
	mfem::Vector center(2);
	center(0) = 0.35;
	center(1) = 0.4;
	pole.omega_p_profile = std::make_shared<GaussianProfile>(3.0, 0.25, center);
	run.model->setLorentzProperties({pole});
	run.open(2);
	expectFullStateAgrees(run, "profiled drude", true);
}

TEST(HesthavenElementLocal, pml_kappa_one_matches_global)
{
	Mesh serial = Mesh::MakeCartesian2D(4, 4, Element::TRIANGLE, true, 1.0, 1.0);
	tagRightHalf(serial, 2);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	mats.emplace(2, Material(1.0, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	PMLProperties pml;
	pml.geom_tags = {2};
	pml.active_axes = {X};
	pml.kappa_max = 1.0;
	pml.alpha_max = 0.0;
	pml.target_reflection = 1e-4;
	pml.grading_order = 2;
	run.model->setPMLProperties({pml});
	run.model->initializePMLProfiles(Mpi::WorldRank(), 2);
	run.open(2);
	expectFullStateAgrees(run, "pml kappa 1", true);
	expectImplicitAgrees(run, "implicit pml kappa 1");
}

TEST(HesthavenElementLocal, pml_kappa_above_one_matches_global)
{
	Mesh serial = Mesh::MakeCartesian2D(4, 4, Element::TRIANGLE, true, 1.0, 1.0);
	tagRightHalf(serial, 2);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	mats.emplace(2, Material(1.0, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	PMLProperties pml;
	pml.geom_tags = {2};
	pml.active_axes = {X, Y};
	pml.kappa_max = 2.0;
	pml.alpha_max = 0.0;
	pml.target_reflection = 1e-4;
	pml.grading_order = 2;
	run.model->setPMLProperties({pml});
	run.model->initializePMLProfiles(Mpi::WorldRank(), 2);
	run.open(2);
	expectFullStateAgrees(run, "pml kappa 2", true);
	expectImplicitAgrees(run, "implicit pml kappa 2");
}

TEST(HesthavenElementLocal, curved_debye_is_not_scaled_twice)
{
	const double eps_inf = 2.0;
	const double sigma = 0.15;
	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	serial.SetCurvature(2);
	GridFunction* nodes = serial.GetNodes();
	ASSERT_NE(nodes, nullptr);
	FiniteElementSpace* node_fes = nodes->FESpace();
	bool moved = false;
	for (int i = 0; i < node_fes->GetNDofs(); ++i) {
		const int ix = node_fes->DofToVDof(i, 0);
		const int iy = node_fes->DofToVDof(i, 1);
		const double x = (*nodes)(ix);
		const double y = (*nodes)(iy);
		if (x > 0.2 && x < 0.8 && y > 0.2 && y < 0.8) {
			(*nodes)(ix) += 0.07;
			moved = true;
			break;
		}
	}
	ASSERT_TRUE(moved);

	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, sigma));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(serial), {}));
	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = 5.0;
	pole.tau_solver = 0.8;
	run.model->setDebyeProperties({pole});
	run.open(2);

	const auto positions = buildDoFPositions(*run.fes);
	const auto lists = initCurvedAndLinearElementsLists(*run.fes, positions);
	ASSERT_GT(lists.first.Size(), 0);
	ASSERT_FALSE(lists.second.empty());
	expectFullStateAgrees(run, "curved debye", true);
}

void expectImplicitAgrees(PairRun& run, const char* name)
{
	const int height = run.global->Height();
	ASSERT_EQ(run.hesthaven->Height(), height) << name;
	Vector x(height);
	for (int i = 0; i < height; ++i) {
		x[i] = 0.1 * std::sin(0.17 * static_cast<double>(i + 1));
	}
	Vector k_global;
	Vector k_hest;
	run.global->ImplicitSolve(0.02, x, k_global);
	run.hesthaven->ImplicitSolve(0.02, x, k_hest);
	const double diff = maxAbsDiff(k_global, k_hest);
	EXPECT_LT(diff, 1e-5) << name << " max |k_global - k_hest| = " << diff;
}

TEST(HesthavenElementLocal, centered_sma_matches_global)
{
	Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	PairRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(smaBoundary(serial), {}));
	run.open(1, 0.0);
	expectFullStateAgrees(run, "centered SMA", false);
}

TEST(HesthavenElementLocal, implicit_conductivity_and_debye_match_global)
{
	{
		Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
		tagRightHalf(serial, 2);
		GeomTagToMaterial mats;
		mats.emplace(1, Material(1.0, 1.0, 0.0));
		mats.emplace(2, Material(2.0, 1.5, 0.2));
		PairRun run;
		run.model = std::make_unique<Model>(
			serial,
			GeomTagToMaterialInfo(mats, {}),
			GeomTagToBoundaryInfo(pecBoundary(serial), {}));
		run.open(1);
		expectImplicitAgrees(run, "implicit conductivity");
	}
	{
		const double eps_inf = 2.0;
		Mesh serial = Mesh::MakeCartesian2D(3, 3, Element::TRIANGLE, true, 1.0, 1.0);
		GeomTagToMaterial mats;
		mats.emplace(1, Material(eps_inf, 1.0, 0.0));
		PairRun run;
		run.model = std::make_unique<Model>(
			serial,
			GeomTagToMaterialInfo(mats, {}),
			GeomTagToBoundaryInfo(pecBoundary(serial), {}));
		DebyeProperties pole;
		pole.geom_tag = 1;
		pole.eps_inf = eps_inf;
		pole.eps_s = 6.0;
		pole.tau_solver = 1.0;
		run.model->setDebyeProperties({pole});
		run.open(1);
		expectImplicitAgrees(run, "implicit debye");
	}
}

TEST(HesthavenElementLocal, interior_boundary_matches_global)
{
	std::ifstream in(maxwellCase("2D_Parser_BdrAndInterior"));
	auto data = nlohmann::json::parse(in);
	data["solver_options"]["order"] = 1;
	auto global_data = data;
	auto hest_data = data;
	global_data["solver_options"]["evolution_operator"] = "global";
	hest_data["solver_options"]["evolution_operator"] = "hesthaven";
	const std::string path = maxwellCase("2D_Parser_BdrAndInterior");
	auto global_solver = driver::buildSolver(global_data, path, true);
	auto hest_solver = driver::buildSolver(hest_data, path, true);
	auto* global = dynamic_cast<GlobalEvolution*>(global_solver.getEvolTDO());
	auto* hest = dynamic_cast<HesthavenEvolution*>(hest_solver.getEvolTDO());
	ASSERT_NE(global, nullptr);
	ASSERT_NE(hest, nullptr);
	ASSERT_EQ(global->Height(), hest->Height());
	Vector x(global->Height());
	for (int i = 0; i < x.Size(); ++i) {
		x[i] = 0.1 * std::sin(0.17 * static_cast<double>(i + 1));
	}
	Vector y_global(x.Size());
	Vector y_hest(x.Size());
	global->Mult(x, y_global);
	hest->Mult(x, y_hest);
	const double diff = maxAbsDiff(y_global, y_hest);
	EXPECT_LT(diff, 1e-8) << "interior boundary max |y_global - y_hest| = " << diff;
}

TEST(HesthavenElementLocal, sgbc_step_matches_global)
{
	std::ifstream in(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
	auto data = nlohmann::json::parse(in);
	data["solver_options"]["order"] = 1;
	data["solver_options"]["time_step"] = 1.0e-3;
	data["solver_options"]["final_time"] = 1.0e-3;
	auto global_data = data;
	auto hest_data = data;
	global_data["solver_options"]["evolution_operator"] = "global";
	hest_data["solver_options"]["evolution_operator"] = "hesthaven";
	const std::string path = maxwellCase("2D_InteriorBoundary_SGBC_Test");
	auto global_solver = driver::buildSolver(global_data, path, true);
	auto hest_solver = driver::buildSolver(hest_data, path, true);
	global_solver.step(false);
	hest_solver.step(false);
	const double diff = maxAbsDiff(
		global_solver.getFields().allDOFs(),
		hest_solver.getFields().allDOFs());
	EXPECT_LT(diff, 1e-6) << "SGBC step max |global - hesthaven| = " << diff;
}

// Regression: promoted SGBC/SMA faces are rank-local. buildGlobalOperator() does
// ExchangeFaceNbrData, so every rank must enter when any rank owns a promoted element.
TEST(HesthavenElementLocal, sgbc_hesthaven_constructs_with_uneven_promotion)
{
	if (Mpi::WorldSize() < 2) {
		GTEST_SKIP() << "Needs >=2 ranks to expose uneven SGBC face ownership.";
	}
	std::ifstream in(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
	auto data = nlohmann::json::parse(in);
	data["solver_options"]["order"] = 1;
	data["solver_options"]["time_step"] = 1.0e-3;
	data["solver_options"]["final_time"] = 1.0e-3;
	data["solver_options"]["upwind_alpha"] = 0.0;
	data["solver_options"]["evolution_operator"] = "hesthaven";
	const std::string path = maxwellCase("2D_InteriorBoundary_SGBC_Test");
	EXPECT_NO_THROW(driver::buildSolver(data, path, true));
}
