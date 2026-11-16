#include <gtest/gtest.h>

#include "SourceFixtures.h"

#include "solver/ProbesManager.h"
#include "solver/SourcesManager.h"
#include <solver/Solver.h>

#include "components/ExporterContextFields.h"
#include "math/PhysicalConstants.h"

#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class ProbesManagerTest : public ::testing::Test {
	
};

TEST_F(ProbesManagerTest, exporterProbe)
{
	Mesh smesh{ Mesh::MakeCartesian1D(20, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes, fields };

	Probes ps;
	ps.exporterProbes = { ExporterProbe{"ProbesManagerTest"} };

	ProbesManager pM{ ps, fes, fields, SolverOptions{} };

	ASSERT_NO_THROW(pM.updateProbes(0.0));
}

TEST_F(ProbesManagerTest, exporterContextFieldsLorentz)
{
	Mesh smesh{ Mesh::MakeCartesian1D(10, 1.0) };
	for (int e = 5; e < smesh.GetNE(); ++e) {
		smesh.SetAttribute(e, 2);
	}
	smesh.SetAttributes();

	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));

	LorentzProperties pole;
	pole.geom_tag = 2;
	pole.eps_inf = 2.0;
	pole.omega_p = 7931762888.007762 / physicalConstants::speedOfLight_SI;
	pole.omega_1 = 1498962290.0 / physicalConstants::speedOfLight_SI;
	pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
	matInfo.lorentz.push_back(pole);

	GeomTagToBoundaryInfo bdrInfo(
		GeomTagToBoundary{{1, BdrCond::SMA}, {2, BdrCond::SMA}},
		GeomTagToInteriorBoundary{});
	Model model(smesh, matInfo, bdrInfo);
	model.setLorentzProperties(matInfo.lorentz);

	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &model.getMesh(), &fec };

	ExporterContextFields ctx(
		model, fes,
		{"eps_r", "mu_r", "sigma_S_per_m", "eps_inf", "fp_Hz",
		 "f1_Hz", "nuc_Hz", "mat_tag"});

	ASSERT_EQ(ctx.fields().size(), 8u);
	EXPECT_FALSE(ctx.isTimeDependent());

	auto find = [&](const std::string& name) -> const mfem::ParGridFunction& {
		for (const auto& f : ctx.fields()) {
			if (f.name == name) { return *f.gf; }
		}
		throw std::runtime_error("missing context field " + name);
	};

	const auto& eps_r = find("eps_r");
	const auto& eps_inf = find("eps_inf");
	const auto& fp = find("fp_Hz");
	const auto& f1 = find("f1_Hz");
	const auto& nuc = find("nuc_Hz");
	const auto& tag = find("mat_tag");

	const double expected_fp = 7931762888.007762 / (2.0 * M_PI);
	const double expected_f1 = 1498962290.0 / (2.0 * M_PI);
	const double expected_nuc = 3185294866.25 / (2.0 * M_PI);

	mfem::Array<int> dofs;
	for (int e = 0; e < fes.GetMesh()->GetNE(); ++e) {
		fes.GetElementDofs(e, dofs);
		const int attr = fes.GetMesh()->GetAttribute(e);
		for (int j = 0; j < dofs.Size(); ++j) {
			const int dof = dofs[j] < 0 ? -1 - dofs[j] : dofs[j];
			EXPECT_NEAR(tag(dof), attr, 1e-12);
			if (attr == 2) {
				EXPECT_NEAR(eps_r(dof), 2.0, 1e-12);
				EXPECT_NEAR(eps_inf(dof), 2.0, 1e-12);
				EXPECT_NEAR(fp(dof), expected_fp, expected_fp * 1e-12);
				EXPECT_NEAR(f1(dof), expected_f1, expected_f1 * 1e-12);
				EXPECT_NEAR(nuc(dof), expected_nuc, expected_nuc * 1e-12);
			} else {
				EXPECT_NEAR(eps_r(dof), 1.0, 1e-12);
				EXPECT_NEAR(eps_inf(dof), 1.0, 1e-12);
				EXPECT_NEAR(fp(dof), 0.0, 1e-12);
				EXPECT_NEAR(f1(dof), 0.0, 1e-12);
				EXPECT_NEAR(nuc(dof), 0.0, 1e-12);
			}
		}
	}
}

TEST_F(ProbesManagerTest, exporterContextFieldsLorentzProfile)
{
	Mesh smesh{ Mesh::MakeCartesian1D(10, 1.0) };
	for (int e = 5; e < smesh.GetNE(); ++e) {
		smesh.SetAttribute(e, 2);
	}
	smesh.SetAttributes();

	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));

	LorentzProperties pole;
	pole.geom_tag = 2;
	pole.eps_inf = 2.0;
	pole.omega_p = 0.0;
	pole.omega_1 = 0.0;
	pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
	mfem::Vector center(1);
	center(0) = 0.5;
	const double peak = 7931762888.007762 / physicalConstants::speedOfLight_SI;
	pole.omega_p_profile =
		std::make_shared<GaussianProfile>(peak, 0.2, center);
	matInfo.lorentz.push_back(pole);

	GeomTagToBoundaryInfo bdrInfo(
		GeomTagToBoundary{{1, BdrCond::SMA}, {2, BdrCond::SMA}},
		GeomTagToInteriorBoundary{});
	Model model(smesh, matInfo, bdrInfo);
	model.setLorentzProperties(matInfo.lorentz);
	const auto* modelPole = model.findLorentz(2);
	ASSERT_NE(modelPole, nullptr);

	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &model.getMesh(), &fec };

	ExporterContextFields ctx(
		model, fes, {"fp_Hz", "nuc_Hz", "eps_inf", "mat_tag"});

	auto find = [&](const std::string& name) -> const mfem::ParGridFunction& {
		for (const auto& f : ctx.fields()) {
			if (f.name == name) { return *f.gf; }
		}
		throw std::runtime_error("missing context field " + name);
	};
	const auto& fp = find("fp_Hz");
	const auto& eps_inf = find("eps_inf");

	double fp_min = std::numeric_limits<double>::max();
	double fp_max = -std::numeric_limits<double>::max();

	mfem::Array<int> dofs;
	mfem::Vector x(1);
	for (int e = 0; e < fes.GetMesh()->GetNE(); ++e) {
		fes.GetElementDofs(e, dofs);
		const int attr = fes.GetMesh()->GetAttribute(e);
		auto* fe = fes.GetFE(e);
		auto* tr = fes.GetElementTransformation(e);
		for (int j = 0; j < dofs.Size(); ++j) {
			const int dof = dofs[j] < 0 ? -1 - dofs[j] : dofs[j];
			const auto& ip = fe->GetNodes().IntPoint(j);
			tr->Transform(ip, x);
			if (attr == 2) {
				const double expected = modelPole->evalOmegaP(x) *
					physicalConstants::speedOfLight_SI / (2.0 * M_PI);
				EXPECT_NEAR(fp(dof), expected,
					std::max(1.0, expected) * 1e-12);
				EXPECT_NEAR(eps_inf(dof), 2.0, 1e-12);
				fp_min = std::min(fp_min, fp(dof));
				fp_max = std::max(fp_max, fp(dof));
			} else {
				EXPECT_NEAR(fp(dof), 0.0, 1e-12);
			}
		}
	}
	EXPECT_LT(fp_min, fp_max);
}

TEST_F(ProbesManagerTest, exporterContextFields2DProfileUniformInY)
{
	Mesh smesh{ Mesh::MakeCartesian2D(
		8, 2, Element::QUADRILATERAL, true, 1.0, 1.0) };
	mfem::Array<int> verts;
	for (int e = 0; e < smesh.GetNE(); ++e) {
		smesh.GetElementVertices(e, verts);
		double xc = 0.0;
		for (int v = 0; v < verts.Size(); ++v) {
			xc += smesh.GetVertex(verts[v])[0];
		}
		xc /= static_cast<double>(verts.Size());
		smesh.SetAttribute(e, xc >= 0.5 ? 2 : 1);
	}
	smesh.SetAttributes();

	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));

	LorentzProperties pole;
	pole.geom_tag = 2;
	pole.eps_inf = 2.0;
	pole.omega_p = 0.0;
	pole.omega_1 = 0.0;
	pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
	mfem::Vector center(1);
	center(0) = 0.5;
	const double peak = 7931762888.007762 / physicalConstants::speedOfLight_SI;
	pole.omega_p_profile =
		std::make_shared<GaussianProfile>(peak, 0.2, center);
	matInfo.lorentz.push_back(pole);

	GeomTagToBoundaryInfo bdrInfo(
		GeomTagToBoundary{{1, BdrCond::SMA}, {2, BdrCond::SMA}},
		GeomTagToInteriorBoundary{});
	Model model(smesh, matInfo, bdrInfo);
	model.setLorentzProperties(matInfo.lorentz);
	const auto* modelPole = model.findLorentz(2);
	ASSERT_NE(modelPole, nullptr);
	EXPECT_NEAR(modelPole->region_min_x, 0.5, 1e-12);
	EXPECT_NEAR(modelPole->region_max_x, 1.0, 1e-12);

	DG_FECollection fec{ 2, 2, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &model.getMesh(), &fec };
	ExporterContextFields ctx(
		model, fes, {"fp_Hz", "nuc_Hz", "eps_inf", "mat_tag"});

	auto find = [&](const std::string& name) -> const mfem::ParGridFunction& {
		for (const auto& f : ctx.fields()) {
			if (f.name == name) { return *f.gf; }
		}
		throw std::runtime_error("missing context field " + name);
	};
	const auto& fp = find("fp_Hz");
	const auto& eps_inf = find("eps_inf");

	std::map<long long, double> fp_by_x;
	mfem::Array<int> dofs;
	mfem::Vector x(2);
	for (int e = 0; e < fes.GetMesh()->GetNE(); ++e) {
		if (fes.GetMesh()->GetAttribute(e) != 2) {
			continue;
		}
		fes.GetElementDofs(e, dofs);
		auto* fe = fes.GetFE(e);
		auto* tr = fes.GetElementTransformation(e);
		for (int j = 0; j < dofs.Size(); ++j) {
			const int dof = dofs[j] < 0 ? -1 - dofs[j] : dofs[j];
			const auto& ip = fe->GetNodes().IntPoint(j);
			tr->Transform(ip, x);
			const double expected = modelPole->evalOmegaP(x) *
				physicalConstants::speedOfLight_SI / (2.0 * M_PI);
			EXPECT_NEAR(fp(dof), expected, std::max(1.0, expected) * 1e-12);
			EXPECT_NEAR(eps_inf(dof), 2.0, 1e-12);

			const long long key = std::llround(x(0) * 1e12);
			const auto it = fp_by_x.find(key);
			if (it != fp_by_x.end()) {
				EXPECT_NEAR(it->second, fp(dof), 1e-12);
			} else {
				fp_by_x[key] = fp(dof);
			}
		}
	}
	EXPECT_GT(fp_by_x.size(), 1u);
}

TEST_F(ProbesManagerTest, exporterContextFieldsApplicability)
{
	auto fieldNames = [](const ExporterContextFields& ctx) {
		std::vector<std::string> names;
		for (const auto& f : ctx.fields()) {
			names.push_back(f.name);
		}
		return names;
	};

	GeomTagToBoundaryInfo bdrInfo(
		GeomTagToBoundary{{1, BdrCond::SMA}, {2, BdrCond::SMA}},
		GeomTagToInteriorBoundary{});

	auto makeMesh = [] {
		Mesh smesh{ Mesh::MakeCartesian1D(10, 1.0) };
		for (int e = 5; e < smesh.GetNE(); ++e) {
			smesh.SetAttribute(e, 2);
		}
		smesh.SetAttributes();
		return smesh;
	};

	{
		auto smesh = makeMesh();
		GeomTagToMaterialInfo matInfo;
		matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
		matInfo.gt2m.emplace(2, Material(1.0, 1.0, 0.0));
		LorentzProperties pole;
		pole.geom_tag = 2;
		pole.eps_inf = 1.0;
		pole.omega_p = 7931762888.007762 / physicalConstants::speedOfLight_SI;
		pole.omega_1 = 0.0;
		pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
		matInfo.lorentz.push_back(pole);
		Model model(smesh, matInfo, bdrInfo);
		model.setLorentzProperties(matInfo.lorentz);
		DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
		ParFiniteElementSpace fes{ &model.getMesh(), &fec };
		ExporterContextFields ctx(model, fes, {});
		EXPECT_EQ(fieldNames(ctx),
			(std::vector<std::string>{
				"eps_r", "mu_r", "mat_tag", "h_m",
				"eps_inf", "fp_Hz", "nuc_Hz"}));
	}

	{
		auto smesh = makeMesh();
		GeomTagToMaterialInfo matInfo;
		matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
		matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));
		LorentzProperties pole;
		pole.geom_tag = 2;
		pole.eps_inf = 2.0;
		pole.omega_p = 7931762888.007762 / physicalConstants::speedOfLight_SI;
		pole.omega_1 = 1498962290.0 / physicalConstants::speedOfLight_SI;
		pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
		matInfo.lorentz.push_back(pole);
		Model model(smesh, matInfo, bdrInfo);
		model.setLorentzProperties(matInfo.lorentz);
		DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
		ParFiniteElementSpace fes{ &model.getMesh(), &fec };
		ExporterContextFields ctx(model, fes, {});
		EXPECT_EQ(fieldNames(ctx),
			(std::vector<std::string>{
				"eps_r", "mu_r", "mat_tag", "h_m",
				"eps_inf", "fp_Hz", "f1_Hz", "nuc_Hz"}));
	}

	{
		auto smesh = makeMesh();
		GeomTagToMaterialInfo matInfo;
		matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
		matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));
		DebyeProperties pole;
		pole.geom_tag = 2;
		pole.eps_inf = 2.0;
		pole.eps_s = 5.0;
		pole.tau_solver = 1e-9 * physicalConstants::speedOfLight_SI;
		matInfo.debye.push_back(pole);
		Model model(smesh, matInfo, bdrInfo);
		model.setDebyeProperties(matInfo.debye);
		DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
		ParFiniteElementSpace fes{ &model.getMesh(), &fec };
		ExporterContextFields ctx(model, fes, {});
		EXPECT_EQ(fieldNames(ctx),
			(std::vector<std::string>{
				"eps_r", "mu_r", "mat_tag", "h_m",
				"eps_inf", "eps_s", "tau_s"}));
	}

	{
		auto smesh = makeMesh();
		GeomTagToMaterialInfo matInfo;
		matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
		matInfo.gt2m.emplace(2, Material(1.0, 1.0, 0.0));
		Model model(smesh, matInfo, bdrInfo);
		PMLProperties pml;
		pml.geom_tags = {2};
		pml.active_axes = {X};
		pml.kappa_max = 1.0;
		pml.alpha_max = 0.0;
		model.setPMLProperties({pml});
		DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
		ParFiniteElementSpace fes{ &model.getMesh(), &fec };
		ExporterContextFields ctx(model, fes, {});
		EXPECT_EQ(fieldNames(ctx),
			(std::vector<std::string>{
				"eps_r", "mu_r", "mat_tag", "h_m",
				"pml_sigma_x", "pml_kappa_x", "pml_alpha_x",
				"pml_depth_x", "pml_mask"}));
	}
}

TEST_F(ProbesManagerTest, exporterContextFieldsSinglePvdFirstFrame)
{
	Mesh smesh{ Mesh::MakeCartesian1D(10, 1.0) };
	for (int e = 5; e < smesh.GetNE(); ++e) {
		smesh.SetAttribute(e, 2);
	}
	smesh.SetAttributes();

	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(2.0, 1.0, 0.0));
	LorentzProperties pole;
	pole.geom_tag = 2;
	pole.eps_inf = 2.0;
	pole.omega_p = 7931762888.007762 / physicalConstants::speedOfLight_SI;
	pole.omega_1 = 1498962290.0 / physicalConstants::speedOfLight_SI;
	pole.gamma = 3185294866.25 / physicalConstants::speedOfLight_SI;
	matInfo.lorentz.push_back(pole);

	GeomTagToBoundaryInfo bdrInfo(
		GeomTagToBoundary{{1, BdrCond::SMA}, {2, BdrCond::SMA}},
		GeomTagToInteriorBoundary{});
	Model model(smesh, matInfo, bdrInfo);
	model.setLorentzProperties(matInfo.lorentz);

	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &model.getMesh(), &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };

	Probes ps;
	ExporterProbe exporter{ "ExporterContextSinglePvdTest" };
	exporter.save_every = 0.1;
	ps.exporterProbes = { exporter };

	SolverOptions opts;
	opts.final_time = 1.0;
	ProbesManager pM{ ps, fes, fields, opts, &model };
	pM.initializeContextFields(model);

	const std::string base = "exports/ParaView/" + getRunModeTag() +
		"/ExporterContextSinglePvdTest";
	std::filesystem::remove_all(base);
	pM.updateProbes(0.0);
	pM.updateProbes(0.1);

	auto readFile = [](const std::string& path) {
		std::ifstream f(path);
		std::stringstream ss;
		ss << f.rdbuf();
		return ss.str();
	};
	const auto cycle0 = readFile(base + "/Cycle000000/proc000000.vtu");
	const auto cycle1 = readFile(base + "/Cycle000001/proc000000.vtu");

	EXPECT_NE(cycle0.find("Name=\"E\""), std::string::npos);
	EXPECT_NE(cycle0.find("Name=\"fp_Hz\""), std::string::npos);
	EXPECT_NE(cycle0.find("Name=\"mat_tag\""), std::string::npos);
	EXPECT_NE(cycle1.find("Name=\"E\""), std::string::npos);
	EXPECT_EQ(cycle1.find("Name=\"fp_Hz\""), std::string::npos);
	EXPECT_EQ(cycle1.find("Name=\"mat_tag\""), std::string::npos);
}

TEST_F(ProbesManagerTest, fieldProbe)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 2, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };
	SourcesManager sM{ buildGaussianInitialField(), fes, fields };

	Probes probes;
	probes.fieldProbes = {
		FieldProbe{E, Z, {0.5}}
	};

	ProbesManager pM{ probes, fes, fields, SolverOptions{} };

	auto tol {1e-2};
	ASSERT_NO_THROW(pM.updateProbes(0.0));
	EXPECT_NO_THROW(pM.getFieldProbe(0));
	for (const auto& [t, f] : pM.getFieldProbe(0).getFieldMovie()){
		EXPECT_NEAR(1.0, f, tol);
	}

}

TEST_F(ProbesManagerTest, morStateProbe_writesTimeAsFirstLine)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };

	MORStateProbe morProbe;
	morProbe.name = "TestMORStateTime";
	morProbe.record_time_start = 2.5;
	morProbe.record_time_final = 2.5;
	morProbe.saves = 1;

	Probes ps;
	ps.morStateProbes = { morProbe };

	SolverOptions opts;
	opts.final_time = 2.5;

	ProbesManager pM{ ps, fes, fields, opts };
	pM.setCaseName("morStateProbeWriteTimeTest");

	// Reproduce the Solver loop pattern: initial call at t=0, then step through with dt=1e-3
	const double dt = 1e-3;
	const double final_time = opts.final_time;
	double time = 0.0;
	pM.updateProbes(time);
	while (time <= final_time - 1e-8 * dt) {
		time = std::min(time + dt, final_time);
		pM.updateProbes(time);
	}

	std::string file_path = "exports/SimulationData/single-core/morStateProbeWriteTimeTest/MORStateProbes/TestMORStateTime/x_0";
	std::ifstream f(file_path);
	ASSERT_TRUE(f.is_open()) << "Expected file at: " << file_path;

	double saved_time;
	f >> saved_time;
	EXPECT_NEAR(2.5, saved_time, 1e-10);
}

TEST_F(ProbesManagerTest, morStateProbe_recordTimesConsecutive)
{
	Mesh smesh{ Mesh::MakeCartesian1D(5, 1.0) };
	ParMesh mesh = ParMesh(MPI_COMM_WORLD, smesh);
	DG_FECollection fec{ 1, 1, BasisType::GaussLobatto };
	ParFiniteElementSpace fes{ &mesh, &fec };
	Fields<ParFiniteElementSpace, ParGridFunction> fields{ fes };

	MORStateProbe morProbe;
	morProbe.name = "TestMORStateConsecutive";
	morProbe.record_time_start = 0.0;
	morProbe.record_time_final = 6.0;
	morProbe.saves = 7;

	Probes ps;
	ps.morStateProbes = { morProbe };

	SolverOptions opts;
	opts.final_time = 6.0;

	ProbesManager pM{ ps, fes, fields, opts };
	pM.setCaseName("morStateProbeConsecutiveTest");

	// Reproduce the Solver loop pattern: initial call at t=0, then step through with dt=1e-3
	const double dt = 1e-3;
	const double final_time = opts.final_time;
	double time = 0.0;
	pM.updateProbes(time);
	while (time <= final_time - 1e-8 * dt) {
		time = std::min(time + dt, final_time);
		pM.updateProbes(time);
	}

	std::string base_path = "exports/SimulationData/single-core/morStateProbeConsecutiveTest/MORStateProbes/TestMORStateConsecutive/";
	for (int i = 0; i < 7; ++i) {
		std::ifstream f(base_path + "x_" + std::to_string(i));
		ASSERT_TRUE(f.is_open()) << "Expected file x_" << i;

		double saved_time;
		f >> saved_time;
		EXPECT_NEAR(static_cast<double>(i), saved_time, dt)
			<< "Expected time " << i << " in file x_" << i;
	}
}
