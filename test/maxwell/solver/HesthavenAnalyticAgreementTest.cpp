#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "TestUtils.h"
#include "components/DGOperatorFactory.h"
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
	const double n = static_cast<double>(verts.Size());
	x /= n;
	y /= n;
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

std::vector<int> elementDofIds(const ParFiniteElementSpace& fes, int elem)
{
	Array<int> dofs;
	fes.GetElementDofs(elem, dofs);
	std::vector<int> ids(static_cast<size_t>(dofs.Size()));
	for (int i = 0; i < dofs.Size(); ++i) {
		ids[static_cast<size_t>(i)] = dofs[i] < 0 ? -1 - dofs[i] : dofs[i];
	}
	return ids;
}

std::vector<int> dofsWithAttribute(const ParFiniteElementSpace& fes, int attribute)
{
	std::vector<int> ids;
	for (int elem = 0; elem < fes.GetMesh()->GetNE(); ++elem) {
		if (fes.GetMesh()->GetAttribute(elem) != attribute) {
			continue;
		}
		const std::vector<int> local = elementDofIds(fes, elem);
		ids.insert(ids.end(), local.begin(), local.end());
	}
	return ids;
}

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

struct DualRun {
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

	void open(int order)
	{
		const int dim = model->getMesh().Dimension();
		fec = std::make_unique<DG_FECollection>(order, dim, BasisType::GaussLobatto);
		fes = std::make_unique<ParFiniteElementSpace>(&model->getMesh(), fec.get());
		const int ndofs = fes->GetNDofs();
		const int n_aux = computePMLAuxSize(model->getPMLProperties(), ndofs, dim)
			+ model->debyeAuxSize(ndofs)
			+ model->lorentzAuxSize(ndofs);
		fields = std::make_unique<Fields<ParFiniteElementSpace, ParGridFunction>>(*fes, n_aux);
		sources_manager = std::make_unique<SourcesManager>(sources, *fes, *fields);
		opts.order = order;
		opts.alpha = 1.0;
		global = std::make_unique<GlobalEvolution>(
			*fes, *model, *sources_manager, opts, probes, 1.0);
		hesthaven = std::make_unique<HesthavenEvolution>(
			*fes, *model, *sources_manager, opts);
	}
};

std::unique_ptr<Model> makeModel(Mesh& mesh, const GeomTagToMaterial& mats)
{
	return std::make_unique<Model>(
		mesh,
		GeomTagToMaterialInfo(mats, {}),
		GeomTagToBoundaryInfo(pecBoundary(mesh), {}));
}

Vector maxwellResidual(DualRun& run, const Vector& state)
{
	const int ndofs = run.fes->GetNDofs();
	if (run.fes->num_face_nbr_dofs != 0) {
		throw std::runtime_error("Analytic residual expects no face-neighbor DOFs.");
	}
	Vector x_eh(6 * ndofs);
	const double* full = state.HostRead();
	for (int i = 0; i < 6 * ndofs; ++i) {
		x_eh[i] = full[i];
	}
	ProblemDescription pd(*run.model, run.probes, run.sources, run.opts);
	DGOperatorFactory<ParFiniteElementSpace> factory(pd, *run.fes);
	auto op = factory.buildGlobalOperator();
	Vector y(6 * ndofs);
	op->Mult(x_eh, y);
	return y;
}

void expectPairAgrees(const Vector& global_y, const Vector& hest_y, const char* name)
{
	ASSERT_EQ(global_y.Size(), hest_y.Size()) << name;
	const double diff = maxAbsDiff(global_y, hest_y);
	EXPECT_LT(diff, 1e-8) << name << " max |global - hesthaven| = " << diff;
}

void thresholdDense(DenseMatrix& mat)
{
	for (int c = 0; c < mat.Width(); ++c) {
		for (int r = 0; r < mat.Height(); ++r) {
			if (std::abs(mat(r, c)) < 1e-8) {
				mat(r, c) = 0.0;
			}
		}
	}
}

DenseMatrix weightedMass(const FiniteElement& fe, ElementTransformation& tr,
                         Coefficient& left, Coefficient& weight, bool drop_tiny)
{
	InverseIntegrator inverse(new MassIntegrator(left));
	DenseMatrix inv_m;
	DenseMatrix mass;
	inverse.AssembleElementMatrix(fe, tr, inv_m);
	MassIntegrator integrator(weight);
	integrator.AssembleElementMatrix(fe, tr, mass);
	DenseMatrix product(inv_m.Height(), mass.Width());
	Mult(inv_m, mass, product);
	if (drop_tiny) {
		thresholdDense(product);
	}
	return product;
}

} // namespace

TEST(HesthavenAnalytic, debye_matches_global_and_single_pole)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps_inf = 2.0;
	const double eps_s = 6.0;
	const double tau = 1.0;
	const double eps_d = eps_s - eps_inf;
	const double e_value = 0.3;
	const double p_value = 0.4;
	const double tol = 1e-9;
	const double e_from_e = -eps_d / (eps_inf * tau);
	const double e_from_p = 1.0 / (eps_inf * tau);

	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = eps_s;
	pole.tau_solver = tau;
	run.model->setDebyeProperties({pole});
	run.open(1);

	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.global->Height(), 9 * ndofs);
	ASSERT_EQ(run.hesthaven->Height(), run.global->Height());

	Vector x0(run.global->Height());
	x0 = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x0[i] = e_value;
	}
	Vector y_global(x0.Size());
	Vector y_hest(x0.Size());
	run.global->Mult(x0, y_global);
	run.hesthaven->Mult(x0, y_hest);
	const Vector y_maxwell = maxwellResidual(run, x0);

	const double* yg = y_global.HostRead();
	const double* yh = y_hest.HostRead();
	const double* ym = y_maxwell.HostRead();
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(yg[i] - ym[i], e_from_e * e_value, tol);
		EXPECT_NEAR(yh[i] - ym[i], e_from_e * e_value, tol);
		EXPECT_NEAR(yg[6 * ndofs + i], (eps_d / tau) * e_value, tol);
		EXPECT_NEAR(yh[6 * ndofs + i], (eps_d / tau) * e_value, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(yg[i], ym[i], tol) << "global index " << i;
		EXPECT_NEAR(yh[i], ym[i], tol) << "hesthaven index " << i;
	}
	EXPECT_LT(maxAbsDiff(y_global, y_hest), tol);

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[6 * ndofs + i] = p_value;
	}
	Vector y1_global(x1.Size());
	Vector y1_hest(x1.Size());
	run.global->Mult(x1, y1_global);
	run.hesthaven->Mult(x1, y1_hest);
	const double* y1g = y1_global.HostRead();
	const double* y1h = y1_hest.HostRead();
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(y1g[i] - yg[i], e_from_p * p_value, tol);
		EXPECT_NEAR(y1h[i] - yh[i], e_from_p * p_value, tol);
		EXPECT_NEAR(y1g[6 * ndofs + i] - yg[6 * ndofs + i], -p_value / tau, tol);
		EXPECT_NEAR(y1h[6 * ndofs + i] - yh[6 * ndofs + i], -p_value / tau, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(y1g[i], yg[i], tol);
		EXPECT_NEAR(y1h[i], yh[i], tol);
	}
	expectPairAgrees(y1_global, y1_hest, "debye with polarization");
}

TEST(HesthavenAnalytic, conductivity_matches_global_and_ohmic_term)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps = 2.0;
	const double sigma = 0.4;
	const double tol = 1e-9;
	Mesh lossy_mesh = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	Mesh ohmic_free_mesh = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial lossy;
	lossy.emplace(1, Material(eps, 1.0, sigma));
	GeomTagToMaterial ohmic_free;
	ohmic_free.emplace(1, Material(eps, 1.0, 0.0));

	DualRun with_sigma;
	with_sigma.model = makeModel(lossy_mesh, lossy);
	with_sigma.open(1);
	DualRun without_sigma;
	without_sigma.model = makeModel(ohmic_free_mesh, ohmic_free);
	without_sigma.open(1);

	const int ndofs = with_sigma.fes->GetNDofs();
	ASSERT_EQ(without_sigma.fes->GetNDofs(), ndofs);
	Vector x(6 * ndofs);
	x = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x[i] = 0.15 * static_cast<double>((i % 4) + 1);
	}

	Vector y_loss_g(x.Size());
	Vector y_loss_h(x.Size());
	Vector y_free_g(x.Size());
	Vector y_free_h(x.Size());
	with_sigma.global->Mult(x, y_loss_g);
	with_sigma.hesthaven->Mult(x, y_loss_h);
	without_sigma.global->Mult(x, y_free_g);
	without_sigma.hesthaven->Mult(x, y_free_h);

	const double* loss_g = y_loss_g.HostRead();
	const double* loss_h = y_loss_h.HostRead();
	const double* free_g = y_free_g.HostRead();
	const double* free_h = y_free_h.HostRead();
	const double* ex = x.HostRead();
	for (int i = 0; i < ndofs; ++i) {
		const double ohmic = -(sigma / eps) * ex[i];
		EXPECT_NEAR(loss_g[i] - free_g[i], ohmic, tol);
		EXPECT_NEAR(loss_h[i] - free_h[i], ohmic, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(loss_g[i], free_g[i], tol) << "global index " << i;
		EXPECT_NEAR(loss_h[i], free_h[i], tol) << "hesthaven index " << i;
	}
	expectPairAgrees(y_loss_g, y_loss_h, "conductive dielectric");
	expectPairAgrees(y_free_g, y_free_h, "same dielectric without conductivity");
}

TEST(HesthavenAnalytic, lorentz_matches_global_and_single_pole)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps_inf = 2.0;
	const double omega_p = 2.0;
	const double omega_1 = 1.0;
	const double gamma = 0.5;
	const double e_value = 0.3;
	const double p_value = 0.2;
	const double j_value = 0.4;
	const double tol = 1e-9;

	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.omega_p = omega_p;
	pole.omega_1 = omega_1;
	pole.gamma = gamma;
	run.model->setLorentzProperties({pole});
	run.open(1);

	const int ndofs = run.fes->GetNDofs();
	const int p_off = 6 * ndofs;
	const int j_off = 9 * ndofs;
	ASSERT_EQ(run.global->Height(), 12 * ndofs);

	Vector x0(run.global->Height());
	x0 = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x0[i] = e_value;
	}
	Vector y0_g(x0.Size());
	Vector y0_h(x0.Size());
	run.global->Mult(x0, y0_g);
	run.hesthaven->Mult(x0, y0_h);
	const Vector y_maxwell = maxwellResidual(run, x0);
	const double* g0 = y0_g.HostRead();
	const double* h0 = y0_h.HostRead();
	const double* ym = y_maxwell.HostRead();
	for (int i = 0; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(g0[i], ym[i], tol) << "global maxwell index " << i;
		EXPECT_NEAR(h0[i], ym[i], tol) << "hesthaven maxwell index " << i;
	}
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(g0[j_off + i], omega_p * omega_p * e_value, tol);
		EXPECT_NEAR(h0[j_off + i], omega_p * omega_p * e_value, tol);
		EXPECT_NEAR(g0[p_off + i], 0.0, tol);
		EXPECT_NEAR(h0[p_off + i], 0.0, tol);
	}

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[p_off + i] = p_value;
		x1[j_off + i] = j_value;
	}
	Vector y1_g(x1.Size());
	Vector y1_h(x1.Size());
	run.global->Mult(x1, y1_g);
	run.hesthaven->Mult(x1, y1_h);
	const double* g1 = y1_g.HostRead();
	const double* h1 = y1_h.HostRead();
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(g1[i] - g0[i], -j_value / eps_inf, tol);
		EXPECT_NEAR(h1[i] - h0[i], -j_value / eps_inf, tol);
		EXPECT_NEAR(g1[p_off + i] - g0[p_off + i], j_value, tol);
		EXPECT_NEAR(h1[p_off + i] - h0[p_off + i], j_value, tol);
		const double j_dot = -omega_1 * omega_1 * p_value - gamma * j_value;
		EXPECT_NEAR(g1[j_off + i] - g0[j_off + i], j_dot, tol);
		EXPECT_NEAR(h1[j_off + i] - h0[j_off + i], j_dot, tol);
	}
	for (int i = ndofs; i < 6 * ndofs; ++i) {
		EXPECT_NEAR(g1[i], g0[i], tol);
		EXPECT_NEAR(h1[i], h0[i], tol);
	}
	expectPairAgrees(y1_g, y1_h, "lorentz pole");
}

TEST(HesthavenAnalytic, drude_has_no_restoring_term)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double omega_p = 1.5;
	const double j_value = 0.25;
	const double p_value = 0.7;
	const double tol = 1e-9;
	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = 1.0;
	pole.omega_p = omega_p;
	pole.omega_1 = 0.0;
	pole.gamma = 0.0;
	run.model->setLorentzProperties({pole});
	run.open(1);

	const int ndofs = run.fes->GetNDofs();
	const int p_off = 6 * ndofs;
	const int j_off = 9 * ndofs;
	Vector x0(run.global->Height());
	x0 = 0.0;
	Vector y0_g(x0.Size());
	Vector y0_h(x0.Size());
	run.global->Mult(x0, y0_g);
	run.hesthaven->Mult(x0, y0_h);

	Vector x1(x0);
	for (int i = 0; i < ndofs; ++i) {
		x1[p_off + i] = p_value;
		x1[j_off + i] = j_value;
	}
	Vector y1_g(x1.Size());
	Vector y1_h(x1.Size());
	run.global->Mult(x1, y1_g);
	run.hesthaven->Mult(x1, y1_h);
	const double* g0 = y0_g.HostRead();
	const double* h0 = y0_h.HostRead();
	const double* g1 = y1_g.HostRead();
	const double* h1 = y1_h.HostRead();
	for (int i = 0; i < ndofs; ++i) {
		EXPECT_NEAR(g1[j_off + i] - g0[j_off + i], 0.0, tol);
		EXPECT_NEAR(h1[j_off + i] - h0[j_off + i], 0.0, tol);
		EXPECT_NEAR(g1[i] - g0[i], -j_value, tol);
		EXPECT_NEAR(h1[i] - h0[i], -j_value, tol);
	}
	expectPairAgrees(y1_g, y1_h, "drude");
}

TEST(HesthavenAnalytic, constant_profile_matches_scalar_on_both_operators)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps_inf = 1.5;
	const double omega_p = 2.2;
	const double omega_1 = 0.8;
	const double gamma = 0.6;
	auto make_scalar = [&]() {
		Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
		GeomTagToMaterial mats;
		mats.emplace(1, Material(eps_inf, 1.0, 0.0));
		DualRun run;
		run.model = makeModel(serial, mats);
		LorentzProperties pole;
		pole.geom_tag = 1;
		pole.eps_inf = eps_inf;
		pole.omega_p = omega_p;
		pole.omega_1 = omega_1;
		pole.gamma = gamma;
		run.model->setLorentzProperties({pole});
		run.open(2);
		return run;
	};
	auto make_profile = [&]() {
		Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
		GeomTagToMaterial mats;
		mats.emplace(1, Material(eps_inf, 1.0, 0.0));
		DualRun run;
		run.model = makeModel(serial, mats);
		LorentzProperties pole;
		pole.geom_tag = 1;
		pole.eps_inf = eps_inf;
		pole.omega_p = omega_p;
		pole.omega_1 = omega_1;
		pole.gamma = gamma;
		pole.omega_p_profile = std::make_shared<ConstantProfile>(omega_p, 2);
		pole.omega_1_profile = std::make_shared<ConstantProfile>(omega_1, 2);
		pole.gamma_profile = std::make_shared<ConstantProfile>(gamma, 2);
		run.model->setLorentzProperties({pole});
		run.open(2);
		return run;
	};

	DualRun scalar = make_scalar();
	DualRun profile = make_profile();
	const int height = scalar.global->Height();
	ASSERT_EQ(profile.global->Height(), height);
	Vector x(height);
	for (int i = 0; i < height; ++i) {
		x[i] = 0.1 * std::sin(0.37 * static_cast<double>(i + 1));
	}
	Vector y_scalar_g(height);
	Vector y_scalar_h(height);
	Vector y_profile_g(height);
	Vector y_profile_h(height);
	scalar.global->Mult(x, y_scalar_g);
	scalar.hesthaven->Mult(x, y_scalar_h);
	profile.global->Mult(x, y_profile_g);
	profile.hesthaven->Mult(x, y_profile_h);

	EXPECT_LT(maxAbsDiff(y_profile_g, y_scalar_g), 1e-10);
	EXPECT_LT(maxAbsDiff(y_profile_h, y_scalar_h), 1e-10);
	expectPairAgrees(y_scalar_g, y_scalar_h, "scalar lorentz");
	expectPairAgrees(y_profile_g, y_profile_h, "constant-profile lorentz");
	expectPairAgrees(y_profile_g, y_scalar_h, "profile global vs scalar hesthaven");
}

TEST(HesthavenAnalytic, gaussian_profile_matches_mass_integrator_on_both_operators)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(1.0, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = 1.0;
	pole.omega_1 = 0.0;
	pole.gamma = 0.35;
	Vector center(2);
	center(0) = 0.25;
	center(1) = 0.30;
	pole.omega_p_profile = std::make_shared<GaussianProfile>(2.4, 0.28, center);
	run.model->setLorentzProperties({pole});
	run.open(2);

	const int ndofs = run.fes->GetNDofs();
	const int j_off = 9 * ndofs;
	const int p_off = 6 * ndofs;
	const auto dofs = elementDofIds(*run.fes, 0);
	const FiniteElement* fe = run.fes->GetFE(0);
	ElementTransformation* tr = run.fes->GetElementTransformation(0);
	const LorentzProperties* stored = run.model->findLorentz(1);
	ASSERT_NE(stored, nullptr);

	ConstantCoefficient one(1.0);
	FunctionCoefficient wp_sq([&](const Vector& x) {
		const double wp = stored->evalOmegaP(x);
		return wp * wp;
	});
	FunctionCoefficient neg_gamma([&](const Vector& x) {
		return -stored->evalGamma(x);
	});
	const DenseMatrix a_wp = weightedMass(*fe, *tr, one, wp_sq, false);
	const DenseMatrix a_nu = weightedMass(*fe, *tr, one, neg_gamma, false);

	Vector x0(run.global->Height());
	x0 = 0.0;
	Vector y0_g(x0.Size());
	Vector y0_h(x0.Size());
	run.global->Mult(x0, y0_g);
	run.hesthaven->Mult(x0, y0_h);

	Vector e_el(static_cast<int>(dofs.size()));
	Vector j_el(static_cast<int>(dofs.size()));
	Vector x1(x0);
	for (int k = 0; k < e_el.Size(); ++k) {
		e_el[k] = 0.2 * static_cast<double>(k + 1);
		j_el[k] = -0.15 * static_cast<double>(k + 2);
		x1[dofs[static_cast<size_t>(k)]] = e_el[k];
		x1[j_off + dofs[static_cast<size_t>(k)]] = j_el[k];
	}
	Vector y1_g(x1.Size());
	Vector y1_h(x1.Size());
	run.global->Mult(x1, y1_g);
	run.hesthaven->Mult(x1, y1_h);

	Vector from_e(e_el.Size());
	Vector from_j(j_el.Size());
	a_wp.Mult(e_el, from_e);
	a_nu.Mult(j_el, from_j);
	const double* g0 = y0_g.HostRead();
	const double* h0 = y0_h.HostRead();
	const double* g1 = y1_g.HostRead();
	const double* h1 = y1_h.HostRead();
	for (int k = 0; k < e_el.Size(); ++k) {
		const int dof = dofs[static_cast<size_t>(k)];
		const double expected = from_e[k] + from_j[k];
		EXPECT_NEAR(g1[j_off + dof] - g0[j_off + dof], expected, 1e-8);
		EXPECT_NEAR(h1[j_off + dof] - h0[j_off + dof], expected, 1e-8);
		EXPECT_NEAR(g1[p_off + dof] - g0[p_off + dof], j_el[k], 1e-10);
		EXPECT_NEAR(h1[p_off + dof] - h0[p_off + dof], j_el[k], 1e-10);
	}
	expectPairAgrees(y1_g, y1_h, "gaussian drude");
}

TEST(HesthavenAnalytic, debye_and_lorentz_regions_match_on_both_operators)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps_inf_d = 2.0;
	const double eps_s_d = 5.0;
	const double tau_d = 0.8;
	const double eps_inf_l = 1.5;
	const double omega_p_l = 2.5;
	const double omega_1_l = 1.2;
	const double gamma_l = 0.4;
	const double tol = 1e-9;

	Mesh serial = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true, 1.0, 1.0);
	tagRightHalf(serial, 2);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf_d, 1.0, 0.0));
	mats.emplace(2, Material(eps_inf_l, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	DebyeProperties debye_pole;
	debye_pole.geom_tag = 1;
	debye_pole.eps_inf = eps_inf_d;
	debye_pole.eps_s = eps_s_d;
	debye_pole.tau_solver = tau_d;
	run.model->setDebyeProperties({debye_pole});
	LorentzProperties lorentz_pole;
	lorentz_pole.geom_tag = 2;
	lorentz_pole.eps_inf = eps_inf_l;
	lorentz_pole.omega_p = omega_p_l;
	lorentz_pole.omega_1 = omega_1_l;
	lorentz_pole.gamma = gamma_l;
	run.model->setLorentzProperties({lorentz_pole});
	run.open(1);

	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.global->Height(), 15 * ndofs);
	const std::vector<int> debye_dofs = dofsWithAttribute(*run.fes, 1);
	const std::vector<int> lorentz_dofs = dofsWithAttribute(*run.fes, 2);
	ASSERT_FALSE(debye_dofs.empty());
	ASSERT_FALSE(lorentz_dofs.empty());

	const double p_debye_val = 0.35;
	const double p_lor_val = 0.25;
	const double j_lor_val = 0.45;
	const int debye_p_off = 6 * ndofs;
	const int lor_p_off = 9 * ndofs;
	const int lor_j_off = 12 * ndofs;

	Vector x0(run.global->Height());
	x0 = 0.0;
	Vector y0_g(x0.Size());
	Vector y0_h(x0.Size());
	run.global->Mult(x0, y0_g);
	run.hesthaven->Mult(x0, y0_h);

	Vector x1(x0);
	for (int dof : debye_dofs) {
		x1[debye_p_off + dof] = p_debye_val;
	}
	for (int dof : lorentz_dofs) {
		x1[lor_p_off + dof] = p_lor_val;
		x1[lor_j_off + dof] = j_lor_val;
	}
	Vector y1_g(x1.Size());
	Vector y1_h(x1.Size());
	run.global->Mult(x1, y1_g);
	run.hesthaven->Mult(x1, y1_h);

	const double* g0 = y0_g.HostRead();
	const double* h0 = y0_h.HostRead();
	const double* g1 = y1_g.HostRead();
	const double* h1 = y1_h.HostRead();
	for (int dof : debye_dofs) {
		EXPECT_NEAR(g1[dof] - g0[dof], p_debye_val / (eps_inf_d * tau_d), tol);
		EXPECT_NEAR(h1[dof] - h0[dof], p_debye_val / (eps_inf_d * tau_d), tol);
		EXPECT_NEAR(g1[debye_p_off + dof] - g0[debye_p_off + dof], -p_debye_val / tau_d, tol);
		EXPECT_NEAR(h1[debye_p_off + dof] - h0[debye_p_off + dof], -p_debye_val / tau_d, tol);
		EXPECT_NEAR(g1[lor_p_off + dof] - g0[lor_p_off + dof], 0.0, tol);
		EXPECT_NEAR(h1[lor_j_off + dof] - h0[lor_j_off + dof], 0.0, tol);
	}
	for (int dof : lorentz_dofs) {
		EXPECT_NEAR(g1[dof] - g0[dof], -j_lor_val / eps_inf_l, tol);
		EXPECT_NEAR(h1[dof] - h0[dof], -j_lor_val / eps_inf_l, tol);
		EXPECT_NEAR(g1[debye_p_off + dof] - g0[debye_p_off + dof], 0.0, tol);
		EXPECT_NEAR(h1[debye_p_off + dof] - h0[debye_p_off + dof], 0.0, tol);
		EXPECT_NEAR(g1[lor_p_off + dof] - g0[lor_p_off + dof], j_lor_val, tol);
		EXPECT_NEAR(h1[lor_p_off + dof] - h0[lor_p_off + dof], j_lor_val, tol);
		const double j_dot = -omega_1_l * omega_1_l * p_lor_val - gamma_l * j_lor_val;
		EXPECT_NEAR(g1[lor_j_off + dof] - g0[lor_j_off + dof], j_dot, tol);
		EXPECT_NEAR(h1[lor_j_off + dof] - h0[lor_j_off + dof], j_dot, tol);
	}
	expectPairAgrees(y1_g, y1_h, "split debye and lorentz");
}

TEST(HesthavenAnalytic, pml_kappa_matches_element_mass_and_global)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	auto run_case = [](double kappa_max, const char* name) {
		Mesh serial = Mesh::MakeCartesian2D(4, 4, Element::TRIANGLE, true, 1.0, 1.0);
		tagRightHalf(serial, 2);
		GeomTagToMaterial mats;
		mats.emplace(1, Material(1.0, 1.0, 0.0));
		mats.emplace(2, Material(1.0, 1.0, 0.0));
		DualRun run;
		run.model = makeModel(serial, mats);
		PMLProperties pml;
		pml.geom_tags = {2};
		pml.active_axes = {X};
		pml.kappa_max = kappa_max;
		pml.alpha_max = 0.0;
		pml.target_reflection = 1e-4;
		pml.grading_order = 2;
		run.model->setPMLProperties({pml});
		run.model->initializePMLProfiles(Mpi::WorldRank(), 2);
		run.open(2);

		int pml_elem = -1;
		for (int elem = 0; elem < run.fes->GetMesh()->GetNE(); ++elem) {
			if (run.fes->GetMesh()->GetAttribute(elem) == 2) {
				pml_elem = elem;
				break;
			}
		}
		ASSERT_GE(pml_elem, 0);
		const auto dofs = elementDofIds(*run.fes, pml_elem);
		const int ndofs = run.fes->GetNDofs();
		const FiniteElement* fe = run.fes->GetFE(pml_elem);
		ElementTransformation* tr = run.fes->GetElementTransformation(pml_elem);
		const PMLProfileData* profiles = run.model->getPMLProfileData();
		ASSERT_NE(profiles, nullptr);

		const bool needs_a = kappa_max > 1.0 + 1e-12;
		ConstantCoefficient one(1.0);
		SCPMLTensorCoefficient c_b(*profiles, X, SCPMLTensorKind::B);
		SCPMLTensorCoefficient c_k(*profiles, X, SCPMLTensorKind::InvKappa);
		SCPMLTensorCoefficient c_a(*profiles, X, SCPMLTensorKind::A);
		DenseMatrix left_b(0);
		DenseMatrix delta(0);
		if (needs_a) {
			left_b = weightedMass(*fe, *tr, c_a, c_b, true);
			DenseMatrix ma_unit = weightedMass(*fe, *tr, c_a, one, true);
			DenseMatrix m_unit = weightedMass(*fe, *tr, one, one, true);
			delta.SetSize(ma_unit.Height(), ma_unit.Width());
			delta = ma_unit;
			delta.Add(-1.0, m_unit);
			thresholdDense(delta);
		} else {
			left_b = weightedMass(*fe, *tr, one, c_b, true);
		}
		const DenseMatrix left_k = weightedMass(*fe, *tr, one, c_k, true);

		Vector x(run.global->Height());
		x = 0.0;
		Vector e_el(static_cast<int>(dofs.size()));
		for (int k = 0; k < e_el.Size(); ++k) {
			e_el[k] = 0.2 * static_cast<double>(k + 1);
			x[dofs[static_cast<size_t>(k)]] = e_el[k];
		}
		Vector y_g(x.Size());
		Vector y_h(x.Size());
		run.global->Mult(x, y_g);
		run.hesthaven->Mult(x, y_h);
		const Vector y_maxwell = maxwellResidual(run, x);

		Vector b_on_e(e_el.Size());
		Vector k_on_e(e_el.Size());
		left_b.Mult(e_el, b_on_e);
		left_k.Mult(e_el, k_on_e);
		Vector curl_local(e_el.Size());
		const double* ym = y_maxwell.HostRead();
		for (int k = 0; k < curl_local.Size(); ++k) {
			curl_local[k] = ym[dofs[static_cast<size_t>(k)]];
		}
		Vector delta_on_curl(e_el.Size());
		delta_on_curl = 0.0;
		if (needs_a) {
			delta.Mult(curl_local, delta_on_curl);
		}

		const double* g = y_g.HostRead();
		const double* h = y_h.HostRead();
		for (int k = 0; k < e_el.Size(); ++k) {
			const int dof = dofs[static_cast<size_t>(k)];
			const double expected_e = -b_on_e[k] + delta_on_curl[k];
			EXPECT_NEAR(g[dof] - ym[dof], expected_e, 1e-7) << name << " global dof " << dof;
			EXPECT_NEAR(h[dof] - ym[dof], expected_e, 1e-7) << name << " hesthaven dof " << dof;
			EXPECT_NEAR(g[6 * ndofs + dof], k_on_e[k], 1e-7) << name << " global P";
			EXPECT_NEAR(h[6 * ndofs + dof], k_on_e[k], 1e-7) << name << " hesthaven P";
		}
		expectPairAgrees(y_g, y_h, name);
	};

	run_case(1.0, "pml kappa 1");
	run_case(2.0, "pml kappa 2");
}

TEST(HesthavenAnalytic, pml_prefix_keeps_debye_tail_on_both_operators)
{
	if (Mpi::WorldSize() > 1) {
		GTEST_SKIP() << "Analytic residual uses an undivided mesh.";
	}

	const double eps_inf = 2.0;
	const double eps_s = 6.0;
	const double tau = 1.0;
	const double p_value = 0.4;
	const double tol = 1e-8;

	Mesh serial = Mesh::MakeCartesian2D(4, 4, Element::TRIANGLE, true, 1.0, 1.0);
	tagRightHalf(serial, 2);
	GeomTagToMaterial mats;
	mats.emplace(1, Material(eps_inf, 1.0, 0.0));
	mats.emplace(2, Material(1.0, 1.0, 0.0));
	DualRun run;
	run.model = makeModel(serial, mats);
	PMLProperties pml;
	pml.geom_tags = {2};
	pml.active_axes = {X};
	pml.kappa_max = 1.0;
	pml.alpha_max = 0.0;
	pml.target_reflection = 1e-4;
	pml.grading_order = 2;
	run.model->setPMLProperties({pml});
	run.model->initializePMLProfiles(Mpi::WorldRank(), 2);
	DebyeProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = eps_inf;
	pole.eps_s = eps_s;
	pole.tau_solver = tau;
	run.model->setDebyeProperties({pole});
	run.open(2);

	const int ndofs = run.fes->GetNDofs();
	ASSERT_EQ(run.global->Height(), 15 * ndofs);
	const std::vector<int> debye_dofs = dofsWithAttribute(*run.fes, 1);
	const std::vector<int> pml_dofs = dofsWithAttribute(*run.fes, 2);
	ASSERT_FALSE(debye_dofs.empty());
	ASSERT_FALSE(pml_dofs.empty());
	const int p_off = 12 * ndofs;

	Vector x0(run.global->Height());
	x0 = 0.0;
	for (int i = 0; i < ndofs; ++i) {
		x0[i] = 0.2;
	}
	Vector y0_g(x0.Size());
	Vector y0_h(x0.Size());
	run.global->Mult(x0, y0_g);
	run.hesthaven->Mult(x0, y0_h);

	Vector x1(x0);
	for (int dof : debye_dofs) {
		x1[p_off + dof] = p_value;
	}
	Vector y1_g(x1.Size());
	Vector y1_h(x1.Size());
	run.global->Mult(x1, y1_g);
	run.hesthaven->Mult(x1, y1_h);

	const double* g0 = y0_g.HostRead();
	const double* h0 = y0_h.HostRead();
	const double* g1 = y1_g.HostRead();
	const double* h1 = y1_h.HostRead();
	for (int dof : pml_dofs) {
		EXPECT_NEAR(g1[dof] - g0[dof], 0.0, tol);
		EXPECT_NEAR(h1[dof] - h0[dof], 0.0, tol);
		EXPECT_NEAR(g1[p_off + dof] - g0[p_off + dof], 0.0, tol);
		EXPECT_NEAR(h1[p_off + dof] - h0[p_off + dof], 0.0, tol);
	}
	for (int dof : debye_dofs) {
		EXPECT_NEAR(g1[dof] - g0[dof], p_value / (eps_inf * tau), tol);
		EXPECT_NEAR(h1[dof] - h0[dof], p_value / (eps_inf * tau), tol);
		EXPECT_NEAR(g1[p_off + dof] - g0[p_off + dof], -p_value / tau, tol);
		EXPECT_NEAR(h1[p_off + dof] - h0[p_off + dof], -p_value / tau, tol);
	}
	expectPairAgrees(y1_g, y1_h, "pml prefix with debye tail");
}
