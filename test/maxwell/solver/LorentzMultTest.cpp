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
			-omega_1 * omega_1 * p_value - 2.0 * gamma * j_value,
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
