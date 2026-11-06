#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "TestUtils.h"
#include "components/DGOperatorFactory.h"
#include "components/LorentzProperties.h"
#include "components/Model.h"
#include "evolution/GlobalEvolution.h"
#include "solver/SourcesManager.h"

using namespace mfem;
using namespace maxwell;

namespace {

GeomTagToBoundary pecEnds(const Mesh& mesh)
{
	GeomTagToBoundary bdr;
	for (int i = 0; i < mesh.GetNBE(); ++i) {
		bdr[mesh.GetBdrAttribute(i)] = BdrCond::PEC;
	}
	return bdr;
}

struct LorentzRun {
	std::unique_ptr<Model> model;
	std::unique_ptr<DG_FECollection> fec;
	std::unique_ptr<ParFiniteElementSpace> fes;
	Sources sources;
	std::unique_ptr<Fields<ParFiniteElementSpace, ParGridFunction>> fields;
	std::unique_ptr<SourcesManager> sources_manager;
	EvolutionOptions opts;
	Probes probes;
	std::unique_ptr<GlobalEvolution> evolution;

	void open(int order)
	{
		fec = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
		fes = std::make_unique<ParFiniteElementSpace>(&model->getMesh(), fec.get());
		fields = std::make_unique<Fields<ParFiniteElementSpace, ParGridFunction>>(
			*fes, model->lorentzAuxSize(fes->GetNDofs()));
		sources_manager = std::make_unique<SourcesManager>(sources, *fes, *fields);
		opts.order = order;
		opts.alpha = 1.0;
		evolution = std::make_unique<GlobalEvolution>(
			*fes, *model, *sources_manager, opts, probes, 1.0);
	}
};

LorentzRun makeLorentzRun(double eps_inf, double omega_p, double omega_1, double gamma)
{
	Mesh serial = Mesh::MakeCartesian1D(2, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	LorentzRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecEnds(serial), {}));
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.omega_p = omega_p;
	pole.omega_1 = omega_1;
	pole.gamma = gamma;
	run.model->setLorentzProperties({pole});
	return run;
}

} // namespace

TEST(LorentzMult, residual_blocks_match_single_pole)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Lorentz residual check uses an undivided mesh.";
	}

	const double eps_inf = 2.0;
	const double omega_p = 2.0;
	const double omega_1 = 1.0;
	const double gamma = 0.5;
	const double e_value = 0.3;
	const double p_value = 0.2;
	const double j_value = 0.4;
	const double tol = 1e-9;

	LorentzRun run = makeLorentzRun(eps_inf, omega_p, omega_1, gamma);
	run.open(1);
	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.fes->num_face_nbr_dofs, 0);
	ASSERT_EQ(run.evolution->Height(), 12 * ndofs);

	const int p_off = 6 * ndofs;
	const int j_off = 9 * ndofs;

	Vector x0(run.evolution->Height());
	x0 = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x0[i] = e_value;
	}
	Vector y0(run.evolution->Height());
	run.evolution->Mult(x0, y0);

	Vector x_eh(6 * ndofs);
	for (int i = 0; i < 6 * ndofs; ++i) {
		x_eh[i] = x0[i];
	}
	ProblemDescription pd(*run.model, run.probes, run.sources, run.opts);
	DGOperatorFactory<ParFiniteElementSpace> factory(pd, *run.fes);
	auto global_op = factory.buildGlobalOperator();
	Vector y_maxwell(6 * ndofs);
	global_op->Mult(x_eh, y_maxwell);

	for (int i = 0; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(y0[i], y_maxwell[i], tol) << "index " << i;
	}
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y0[p_off + i], 0.0, tol);
		EXPECT_NEAR(y0[j_off + i], omega_p * omega_p * e_value, tol);
	}

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[p_off + i] = p_value;
		x1[j_off + i] = j_value;
	}
	Vector y1(run.evolution->Height());
	run.evolution->Mult(x1, y1);

	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1[i] - y0[i], -j_value / eps_inf, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(y1[i] - y0[i], 0.0, tol);
	}
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1[p_off + i] - y0[p_off + i], j_value, tol);
		EXPECT_NEAR(
			y1[j_off + i] - y0[j_off + i],
			-omega_1 * omega_1 * p_value - gamma * j_value,
			tol);
	}
}

TEST(LorentzMult, cold_plasma_has_no_restoring_term)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Lorentz residual check uses an undivided mesh.";
	}

	const double omega_p = 1.5;
	const double j_value = 0.25;
	const double p_value = 0.7;
	const double tol = 1e-9;
	LorentzRun run = makeLorentzRun(1.0, omega_p, 0.0, 0.0);
	run.open(1);
	const int ndofs = run.fes->GetNDofs();
	const int p_off = 6 * ndofs;
	const int j_off = 9 * ndofs;

	Vector x0(run.evolution->Height());
	x0 = 0.0;
	Vector y0(run.evolution->Height());
	run.evolution->Mult(x0, y0);

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[p_off + i] = p_value;
		x1[j_off + i] = j_value;
	}
	Vector y1(run.evolution->Height());
	run.evolution->Mult(x1, y1);

	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1[j_off + i] - y0[j_off + i], 0.0, tol);
		EXPECT_NEAR(y1[i] - y0[i], -j_value, tol);
	}
}

TEST(LorentzMult, constant_profile_mass_integrator_matches_scalar_pole)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Lorentz residual check uses an undivided mesh.";
	}

	const double eps_inf = 1.5;
	const double omega_p = 2.2;
	const double omega_1 = 0.8;
	const double gamma = 0.6;

	LorentzRun scalar_run = makeLorentzRun(eps_inf, omega_p, omega_1, gamma);
	scalar_run.open(2);

	LorentzRun profile_run = makeLorentzRun(eps_inf, omega_p, omega_1, gamma);
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.omega_p = omega_p;
	pole.omega_1 = omega_1;
	pole.gamma = gamma;
	pole.omega_p_profile = std::make_shared<ConstantProfile>(omega_p);
	pole.omega_1_profile = std::make_shared<ConstantProfile>(omega_1);
	pole.gamma_profile = std::make_shared<ConstantProfile>(gamma);
	profile_run.model->setLorentzProperties({pole});
	profile_run.open(2);

	const int h = scalar_run.evolution->Height();
	Vector x(h), y_scalar(h), y_profile(h);
	for (int i = 0; i < h; ++i) {
		x[i] = 0.1 * std::sin(0.37 * static_cast<double>(i + 1));
	}
	scalar_run.evolution->Mult(x, y_scalar);
	profile_run.evolution->Mult(x, y_profile);

	for (int i = 0; i < h; ++i) {
		EXPECT_NEAR(y_profile[i], y_scalar[i], 1e-11) << "index " << i;
	}
}

TEST(LorentzMult, ramc61km_stagnation_profile_matches_mass_integrator_action)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Lorentz residual check uses an undivided mesh.";
	}

	// Build 2-element mesh on [0, 1.0] where el 0 is [0, 0.988] (vacuum, tag 1)
	// and el 1 is [0.988, 1.0] (stagnation-line Drude plasma, tag 2).
	Mesh serial(1, 3, 2, 2);
	const double vx0 = 0.0, vx1 = 0.988, vx2 = 1.0;
	serial.AddVertex(&vx0);
	serial.AddVertex(&vx1);
	serial.AddVertex(&vx2);
	serial.AddSegment(0, 1, 1);
	serial.AddSegment(1, 2, 2);
	Element* b0 = serial.NewElement(Geometry::POINT);
	int v0 = 0;
	b0->SetVertices(&v0);
	b0->SetAttribute(1);
	serial.AddBdrElement(b0);
	Element* b1 = serial.NewElement(Geometry::POINT);
	int v2 = 2;
	b1->SetVertices(&v2);
	b1->SetAttribute(2);
	serial.AddBdrElement(b1);
	serial.FinalizeTopology();
	serial.Finalize();

	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	mats.emplace(2, Material(1.0, 1.0, 0.0));

	LorentzRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecEnds(serial), {}));

	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	LorentzProperties pole;
	pole.geom_tag = 2;
	pole.eps_inf = 1.0;
	pole.omega_p = 1.0;
	pole.omega_1 = 0.0;
	pole.gamma = 1.0;
	pole.omega_p_profile = std::make_shared<Ramc61kmStagnationProfile>(
		Ramc61kmPreset::OmegaP, inv_c, 1.0, 0, 1.0);
	pole.gamma_profile = std::make_shared<Ramc61kmStagnationProfile>(
		Ramc61kmPreset::NuC_Gupta, inv_c, 1.0, 0, 1.0);
	run.model->setLorentzProperties({pole});
	run.open(2);

	const int ndofs = run.fes->GetNDofs();
	const int p_off = 6 * ndofs;
	const int j_off = 9 * ndofs;

	// Assemble reference element matrix on el 1 using MFEM MassIntegrator + FunctionCoefficient directly.
	ConstantCoefficient one_coeff(1.0);
	InverseIntegrator m_inv_integ(new MassIntegrator(one_coeff));
	FunctionCoefficient wp_sq_fc([&](const Vector& p) {
		Vector s(1);
		s(0) = p(0) - 1.0;
		const double wp = ramc61kmStagnationOmegaP(s) * inv_c;
		return wp * wp;
	});
	FunctionCoefficient nu_c_fc([&](const Vector& p) {
		Vector s(1);
		s(0) = p(0) - 1.0;
		return -ramc61kmStagnationNuC_Gupta(s) * inv_c;
	});
	MassIntegrator mi_wp(wp_sq_fc);
	MassIntegrator mi_nu(nu_c_fc);

	const FiniteElement* fe = run.fes->GetFE(1);
	ElementTransformation* tr = run.fes->GetElementTransformation(1);
	DenseMatrix el_m_inv, el_m_wp, el_m_nu, a_wp(fe->GetDof()), a_nu(fe->GetDof());
	m_inv_integ.AssembleElementMatrix(*fe, *tr, el_m_inv);
	mi_wp.AssembleElementMatrix(*fe, *tr, el_m_wp);
	mi_nu.AssembleElementMatrix(*fe, *tr, el_m_nu);
	Mult(el_m_inv, el_m_wp, a_wp);
	Mult(el_m_inv, el_m_nu, a_nu);

	Array<int> el1_dofs;
	run.fes->GetElementDofs(1, el1_dofs);

	Vector x0(run.evolution->Height());
	x0 = 0.0;
	Vector y0(run.evolution->Height());
	run.evolution->Mult(x0, y0);

	// Set non-uniform E_x and J_x on element 1 and check dJ_x/dt block action.
	Vector x1(x0);
	Vector e_el(el1_dofs.Size()), j_el(el1_dofs.Size());
	for (int k = 0; k < el1_dofs.Size(); ++k) {
		e_el[k] = 0.4 * (k + 1);
		j_el[k] = -0.25 * (k + 2);
		x1[el1_dofs[k]] = e_el[k];
		x1[j_off + el1_dofs[k]] = j_el[k];
	}
	Vector y1(run.evolution->Height());
	run.evolution->Mult(x1, y1);

	Vector dj_expected(el1_dofs.Size());
	a_wp.Mult(e_el, dj_expected);
	a_nu.AddMult(j_el, dj_expected);

	for (int k = 0; k < el1_dofs.Size(); ++k) {
		EXPECT_NEAR(y1[j_off + el1_dofs[k]] - y0[j_off + el1_dofs[k]], dj_expected[k], 1e-8);
		EXPECT_NEAR(y1[p_off + el1_dofs[k]] - y0[p_off + el1_dofs[k]], j_el[k], 1e-11);
	}
	// Element 0 (vacuum) has zero J update.
	Array<int> el0_dofs;
	run.fes->GetElementDofs(0, el0_dofs);
	for (int k = 0; k < el0_dofs.Size(); ++k) {
		EXPECT_NEAR(y1[j_off + el0_dofs[k]], 0.0, 1e-12);
	}
}

TEST(LorentzMult, profile_stretches_across_arbitrary_1d_region_lengths)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Lorentz residual check uses an undivided mesh.";
	}

	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;

	// Test stretching the 12 mm RAM-C profile across a 24 mm region [0.976, 1.000]
	// and a 2.0 m region [0.5, 2.5].
	for (const auto& [x_start, x_end] : std::vector<std::pair<double, double>>{{0.976, 1.000}, {0.5, 2.5}}) {
		Mesh serial(1, 3, 2, 2);
		const double vx0 = 0.0, vx1 = x_start, vx2 = x_end;
		serial.AddVertex(&vx0);
		serial.AddVertex(&vx1);
		serial.AddVertex(&vx2);
		serial.AddSegment(0, 1, 1);
		serial.AddSegment(1, 2, 2);
		Element* b0 = serial.NewElement(Geometry::POINT);
		int v0 = 0;
		b0->SetVertices(&v0);
		b0->SetAttribute(1);
		serial.AddBdrElement(b0);
		Element* b1 = serial.NewElement(Geometry::POINT);
		int v2 = 2;
		b1->SetVertices(&v2);
		b1->SetAttribute(2);
		serial.AddBdrElement(b1);
		serial.FinalizeTopology();
		serial.Finalize();

		GeomTagToMaterial mats;
		mats.emplace(1, Material(1.0, 1.0, 0.0));
		mats.emplace(2, Material(1.0, 1.0, 0.0));

		Model model(
			serial,
			GeomTagToMaterialInfo(mats, {}),
			GeomTagToBoundaryInfo(pecEnds(serial), {}));

		LorentzProperties pole;
		pole.geom_tag = 2;
		pole.eps_inf = 1.0;
		pole.omega_p = 1.0;
		pole.gamma = 1.0;
		pole.omega_p_profile = std::make_shared<Ramc61kmStagnationProfile>(
			Ramc61kmPreset::OmegaP, inv_c);
		pole.gamma_profile = std::make_shared<Ramc61kmStagnationProfile>(
			Ramc61kmPreset::NuC_Gupta, inv_c);
		model.setLorentzProperties({pole});

		const LorentzProperties* stored = model.findLorentz(2);
		ASSERT_NE(stored, nullptr);
		EXPECT_DOUBLE_EQ(stored->region_min_x, x_start);
		EXPECT_DOUBLE_EQ(stored->region_max_x, x_end);

		for (double xi : {0.0, 0.25, 0.5, 0.75, 1.0}) {
			Vector pt_mesh(1);
			pt_mesh(0) = x_start + xi * (x_end - x_start);
			Vector pt_native(1);
			pt_native(0) = -0.012 * (1.0 - xi);
			EXPECT_NEAR(stored->evalOmegaP(pt_mesh), ramc61kmStagnationOmegaP(pt_native) * inv_c, 1e-12);
			EXPECT_NEAR(stored->evalGamma(pt_mesh), ramc61kmStagnationNuC_Gupta(pt_native) * inv_c, 1e-12);
		}
	}
}
