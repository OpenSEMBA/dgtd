#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "TestUtils.h"
#include "components/DGOperatorFactory.h"
#include "components/DebyeProperties.h"
#include "components/Model.h"
#include "components/PMLProperties.h"
#include "evolution/GlobalEvolution.h"
#include "solver/SourcesManager.h"

using namespace mfem;
using namespace maxwell;

namespace {

double hostMaxAbsDiff(const Vector& a, const Vector& b, int begin, int end)
{
	const double* ah = a.HostRead();
	const double* bh = b.HostRead();
	double m = 0.0;
	for (int i = begin; i < end; ++i) {
		m = std::max(m, std::abs(ah[i] - bh[i]));
	}
	return m;
}

double hostMaxAbs(const Vector& a, int begin, int end)
{
	const double* ah = a.HostRead();
	double m = 0.0;
	for (int i = begin; i < end; ++i) {
		m = std::max(m, std::abs(ah[i]));
	}
	return m;
}

GeomTagToBoundary pecEnds(const Mesh& mesh)
{
	GeomTagToBoundary bdr;
	for (int i = 0; i < mesh.GetNBE(); ++i) {
		bdr[mesh.GetBdrAttribute(i)] = BdrCond::PEC;
	}
	return bdr;
}

struct DebyeRun {
	std::unique_ptr<Model> model;
	std::unique_ptr<DG_FECollection> fec;
	std::unique_ptr<ParFiniteElementSpace> fes;
	Sources sources;
	std::unique_ptr<Fields<ParFiniteElementSpace, ParGridFunction>> fields;
	std::unique_ptr<SourcesManager> sources_manager;
	EvolutionOptions opts;
	Probes probes;
	std::unique_ptr<GlobalEvolution> evolution;

	void open(int order, int dim)
	{
		fec = std::make_unique<DG_FECollection>(order, dim, BasisType::GaussLobatto);
		fes = std::make_unique<ParFiniteElementSpace>(&model->getMesh(), fec.get());
		const int n_aux = computePMLAuxSize(
			model->getPMLProperties(), fes->GetNDofs(), dim);
		fields = std::make_unique<Fields<ParFiniteElementSpace, ParGridFunction>>(
			*fes, n_aux + model->debyeAuxSize(fes->GetNDofs()));
		sources_manager = std::make_unique<SourcesManager>(sources, *fes, *fields);
		opts.order = order;
		opts.alpha = 1.0;
		evolution = std::make_unique<GlobalEvolution>(
			*fes, *model, *sources_manager, opts, probes, 1.0);
	}
};

std::vector<int> dofsWithAttribute(const ParFiniteElementSpace& fes, int attribute)
{
	std::vector<int> ids;
	mfem::Array<int> dofs;
	for (int el = 0; el < fes.GetMesh()->GetNE(); ++el) {
		if (fes.GetMesh()->GetAttribute(el) != attribute) {
			continue;
		}
		fes.GetElementDofs(el, dofs);
		for (int j = 0; j < dofs.Size(); ++j) {
			int dof = dofs[j];
			if (dof < 0) {
				dof = -1 - dof;
			}
			ids.push_back(dof);
		}
	}
	return ids;
}

} // namespace

TEST(DebyeMult, residual_blocks_match_single_pole)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Debye residual check uses an undivided mesh.";
	}

	const int order = 1;
	const double eps_inf = 2.0;
	const double eps_s = 6.0;
	const double tau = 1.0;
	const double eps_d = eps_s - eps_inf;
	const double e_value = 0.3;
	const double p_value = 0.4;
	const double tol = 1e-9;

	Mesh serial = Mesh::MakeCartesian1D(2, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	DebyeRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecEnds(serial), {}));
	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = eps_s;
	pole.tau_solver = tau;
	run.model->setDebyeProperties({pole});
	run.open(order, 1);
	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.fes->num_face_nbr_dofs, 0);
	ASSERT_EQ(run.evolution->Height(), 9 * ndofs);

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

	const double e_from_e = -eps_d / (eps_inf * tau);
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y0[i] - y_maxwell[i], e_from_e * e_value, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(y0[i], y_maxwell[i], tol) << "index " << i;
	}
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y0[6 * ndofs + i], (eps_d / tau) * e_value, tol);
	}
	EXPECT_LT(hostMaxAbs(y0, 6 * ndofs + ndofs, 9 * ndofs), tol);

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[6 * ndofs + i] = p_value;
	}
	Vector y1(run.evolution->Height());
	run.evolution->Mult(x1, y1);

	const double e_from_p = 1.0 / (eps_inf * tau);
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1[i] - y0[i], e_from_p * p_value, tol);
	}
	EXPECT_LT(hostMaxAbsDiff(y1, y0, ndofs, 6 * ndofs), tol);
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1[6 * ndofs + i] - y0[6 * ndofs + i], -p_value / tau, tol);
	}
	EXPECT_LT(hostMaxAbsDiff(y1, y0, 7 * ndofs, 9 * ndofs), tol);
}

TEST(DebyeMult, pml_prefix_keeps_debye_in_the_tail)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Debye residual check uses an undivided mesh.";
	}

	const int order = 1;
	const double eps_inf = 2.0;
	const double eps_s = 6.0;
	const double tau = 1.0;
	const double p_value = 0.4;
	const double tol = 1e-9;

	Mesh serial = Mesh::MakeCartesian1D(4, 1.0);
	serial.SetAttribute(0, 2);
	serial.SetAttribute(1, 2);
	serial.SetAttributes();

	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	mats.emplace(2, Material(1.0, 1.0, 0.0));
	DebyeRun run;
	run.model = std::make_unique<Model>(
		serial,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecEnds(serial), {}));

	PMLProperties pml;
	pml.geom_tags = {2};
	pml.active_axes = {X};
	pml.kappa_max = 1.0;
	pml.alpha_max = 0.0;
	pml.target_reflection = 1e-4;
	pml.grading_order = 1;
	run.model->setPMLProperties({pml});
	run.model->initializePMLProfiles(Mpi::WorldRank(), order);

	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = eps_s;
	pole.tau_solver = tau;
	run.model->setDebyeProperties({pole});
	run.open(order, 1);

	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.evolution->Height(), 15 * ndofs);
	const std::vector<int> debye_dofs = dofsWithAttribute(*run.fes, 1);
	const std::vector<int> pml_dofs = dofsWithAttribute(*run.fes, 2);
	ASSERT_FALSE(debye_dofs.empty());
	ASSERT_FALSE(pml_dofs.empty());

	Vector x0(run.evolution->Height());
	x0 = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x0[i] = 0.2;
	}
	Vector y0(run.evolution->Height());
	run.evolution->Mult(x0, y0);

	Vector x1(x0);
	const int p_off = 12 * ndofs;
	for (int dof : debye_dofs) {
		x1[p_off + dof] = p_value;
	}
	Vector y1(run.evolution->Height());
	run.evolution->Mult(x1, y1);

	EXPECT_LT(hostMaxAbsDiff(y1, y0, 6 * ndofs, 12 * ndofs), tol);
	EXPECT_LT(hostMaxAbsDiff(y1, y0, ndofs, 6 * ndofs), tol);
	for (int dof : pml_dofs) {
		EXPECT_NEAR(y1[dof] - y0[dof], 0.0, tol);
		EXPECT_NEAR(y1[p_off + dof] - y0[p_off + dof], 0.0, tol);
	}
	for (int dof : debye_dofs) {
		EXPECT_NEAR(y1[dof] - y0[dof], p_value / (eps_inf * tau), tol);
		EXPECT_NEAR(y1[p_off + dof] - y0[p_off + dof], -p_value / tau, tol);
	}
	EXPECT_LT(hostMaxAbsDiff(y1, y0, p_off + ndofs, 15 * ndofs), tol);
}
