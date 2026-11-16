#include "evolution/HesthavenAuxiliary.h"

#include "components/DGOperatorFactory.h"
#include "components/SCPMLLayout.h"
#include "math/EigenMfemTools.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef SEMBA_DGTD_ENABLE_OPENMP
#include <omp.h>
#endif

namespace maxwell {

namespace {

const HesthavenLinearMaterial kVacuumMaterial{};

std::vector<int> elementDofs(mfem::FiniteElementSpace& fes, int elem)
{
	mfem::Array<int> vdofs;
	fes.GetElementDofs(elem, vdofs);
	std::vector<int> dofs(static_cast<size_t>(vdofs.Size()));
	for (int i = 0; i < vdofs.Size(); ++i) {
		const int dof = vdofs[i];
		dofs[static_cast<size_t>(i)] = dof < 0 ? -1 - dof : dof;
	}
	return dofs;
}

DynamicMatrix toEigenThreshold(const mfem::DenseMatrix& mat, bool threshold)
{
	DynamicMatrix out(mat.Height(), mat.Width());
	for (int c = 0; c < mat.Width(); ++c) {
		for (int r = 0; r < mat.Height(); ++r) {
			double value = mat(r, c);
			if (threshold && std::abs(value) < 1e-8) {
				value = 0.0;
			}
			out(r, c) = value;
		}
	}
	return out;
}

void assembleMass(const mfem::FiniteElement& fe, mfem::ElementTransformation& tr,
                  mfem::Coefficient& coeff, mfem::DenseMatrix& mass)
{
	mfem::MassIntegrator integrator(coeff);
	integrator.AssembleElementMatrix(fe, tr, mass);
}

void assembleInverseMass(const mfem::FiniteElement& fe, mfem::ElementTransformation& tr,
                         mfem::Coefficient& coeff, mfem::DenseMatrix& inverse)
{
	mfem::InverseIntegrator integrator(new mfem::MassIntegrator(coeff));
	integrator.AssembleElementMatrix(fe, tr, inverse);
}

DynamicMatrix weighted(const mfem::DenseMatrix& inverse, mfem::Coefficient& coeff,
                       const mfem::FiniteElement& fe, mfem::ElementTransformation& tr,
                       bool threshold)
{
	mfem::DenseMatrix mass;
	mfem::DenseMatrix product;
	assembleMass(fe, tr, coeff, mass);
	product.SetSize(inverse.Height(), mass.Width());
	mfem::Mult(inverse, mass, product);
	return toEigenThreshold(product, threshold);
}

void addMatvec(const DynamicMatrix& matrix, double scale,
               const double* in, int in_offset,
               double* out, int out_offset,
               const std::vector<int>& dofs, int /*ndofs*/)
{
	if (matrix.size() == 0 || scale == 0.0) {
		return;
	}
	const int n = static_cast<int>(dofs.size());
#ifdef SEMBA_DGTD_ENABLE_OPENMP
	// One scratch pair per OpenMP thread — avoids per-call heap traffic.
#endif
	thread_local Eigen::VectorXd x;
	thread_local Eigen::VectorXd y;
	x.resize(n);
	y.resize(n);
	for (int i = 0; i < n; ++i) {
		x[i] = in[in_offset + dofs[static_cast<size_t>(i)]];
	}
	y.noalias() = matrix * x;
	if (scale != 1.0) {
		y *= scale;
	}
	for (int i = 0; i < n; ++i) {
		out[out_offset + dofs[static_cast<size_t>(i)]] += y[i];
	}
}

void accumulateDenseBlock(mfem::SparseMatrix& A, const DynamicMatrix& M, double scale,
                          const std::vector<int>& dofs, int row_off, int col_off)
{
	if (M.size() == 0 || scale == 0.0) {
		return;
	}
	const int n = static_cast<int>(dofs.size());
	for (int j = 0; j < n; ++j) {
		const int gj = col_off + dofs[static_cast<size_t>(j)];
		for (int i = 0; i < n; ++i) {
			const double v = scale * M(i, j);
			if (v != 0.0) {
				A.Add(row_off + dofs[static_cast<size_t>(i)], gj, v);
			}
		}
	}
}

HesthavenPMLElement buildPmlElement(mfem::ParFiniteElementSpace& fes, Model& model,
                                    int elem, bool needs_a, const PMLProfileData& profiles)
{
	HesthavenPMLElement term;
	term.dofs = elementDofs(fes, elem);
	term.needs_a = needs_a;

	const auto* fe = fes.GetFE(elem);
	mfem::ElementTransformation* tr = fes.GetElementTransformation(elem);
	const int attr = fes.GetMesh()->GetAttribute(elem);
	const auto mat_it = model.getGeomTagToMaterial().find(attr);
	const double eps = mat_it == model.getGeomTagToMaterial().end() ? 1.0 : mat_it->second.getPermittivity();
	const double mu = mat_it == model.getGeomTagToMaterial().end() ? 1.0 : mat_it->second.getPermeability();

	mfem::ConstantCoefficient eps_coeff(eps);
	mfem::ConstantCoefficient mu_coeff(mu);
	mfem::ConstantCoefficient one(1.0);
	mfem::DenseMatrix inv_e;
	mfem::DenseMatrix inv_h;
	mfem::DenseMatrix mass_one;
	assembleInverseMass(*fe, *tr, eps_coeff, inv_e);
	assembleInverseMass(*fe, *tr, mu_coeff, inv_h);
	assembleMass(*fe, *tr, one, mass_one);

	for (int axis = 0; axis < 3; ++axis) {
		const Direction dir = static_cast<Direction>(axis);
		SCPMLTensorCoefficient c_a(profiles, dir, SCPMLTensorKind::A);
		SCPMLTensorCoefficient c_b(profiles, dir, SCPMLTensorKind::B);
		SCPMLTensorCoefficient c_c(profiles, dir, SCPMLTensorKind::C);
		SCPMLTensorCoefficient c_d(profiles, dir, SCPMLTensorKind::D);
		SCPMLTensorCoefficient c_invk(profiles, dir, SCPMLTensorKind::InvKappa);

		mfem::DenseMatrix left_e = inv_e;
		mfem::DenseMatrix left_h = inv_h;
		if (needs_a) {
			mfem::DenseMatrix inv_a;
			assembleInverseMass(*fe, *tr, c_a, inv_a);
			left_e = inv_a;
			left_h = inv_a;
			mfem::DenseMatrix scaled(inv_a.Height(), mass_one.Width());
			mfem::DenseMatrix unit_e(inv_e.Height(), mass_one.Width());
			mfem::Mult(inv_a, mass_one, scaled);
			mfem::Mult(inv_e, mass_one, unit_e);
			scaled.Add(-1.0, unit_e);
			term.delta[static_cast<size_t>(axis)] = toEigenThreshold(scaled, true);
		}

		term.blocks[static_cast<size_t>(axis)][0] = weighted(left_e, c_b, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][1] = weighted(left_h, c_b, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][2] = weighted(left_e, c_c, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][3] = weighted(left_h, c_c, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][4] = weighted(inv_e, c_d, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][5] = weighted(inv_h, c_d, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][6] = weighted(inv_e, c_invk, *fe, *tr, true);
		term.blocks[static_cast<size_t>(axis)][7] = weighted(inv_h, c_invk, *fe, *tr, true);
	}
	return term;
}

HesthavenLorentzTerm buildLorentzTerm(mfem::ParFiniteElementSpace& fes, int elem,
                                      const LorentzProperties& pole)
{
	HesthavenLorentzTerm term;
	term.dofs = elementDofs(fes, elem);
	term.e_from_j = -1.0 / pole.eps_inf;
	term.spatial = pole.hasSpatialProfile();
	if (!term.spatial) {
		term.j_from_e = pole.omega_p * pole.omega_p;
		term.j_from_p = -pole.omega_1 * pole.omega_1;
		term.j_from_j = -pole.gamma;
		return term;
	}

	const auto* fe = fes.GetFE(elem);
	mfem::ElementTransformation* tr = fes.GetElementTransformation(elem);
	mfem::ConstantCoefficient one(1.0);
	mfem::DenseMatrix inv_m;
	assembleInverseMass(*fe, *tr, one, inv_m);

	mfem::FunctionCoefficient wp_sq([&pole](const mfem::Vector& x) {
		const double wp = pole.evalOmegaP(x);
		return wp * wp;
	});
	term.j_from_e_block = weighted(inv_m, wp_sq, *fe, *tr, false);

	term.has_p_block = pole.omega_1_profile != nullptr || pole.omega_1 != 0.0;
	if (term.has_p_block) {
		mfem::FunctionCoefficient neg_w1([&pole](const mfem::Vector& x) {
			const double w1 = pole.evalOmega1(x);
			return -w1 * w1;
		});
		term.j_from_p_block = weighted(inv_m, neg_w1, *fe, *tr, false);
	}

	mfem::FunctionCoefficient neg_gamma([&pole](const mfem::Vector& x) {
		return -pole.evalGamma(x);
	});
	term.j_from_j_block = weighted(inv_m, neg_gamma, *fe, *tr, false);
	return term;
}

void applyPmlElement(const HesthavenPMLElement& term, int ndofs,
                     const double* in, double* out)
{
	for (int axis = 0; axis < 3; ++axis) {
		if (!term.needs_a) {
			continue;
		}
		const int e_off = axis * ndofs;
		const int h_off = (3 + axis) * ndofs;
		addMatvec(term.delta[static_cast<size_t>(axis)], 1.0,
		          out, e_off, out, e_off, term.dofs, ndofs);
		addMatvec(term.delta[static_cast<size_t>(axis)], 1.0,
		          out, h_off, out, h_off, term.dofs, ndofs);
	}
	for (int axis = 0; axis < 3; ++axis) {
		const int e_off = axis * ndofs;
		const int h_off = (3 + axis) * ndofs;
		const int pe = 6 * ndofs + axis * ndofs;
		const int ph = 6 * ndofs + (3 + axis) * ndofs;
		const auto& block = term.blocks[static_cast<size_t>(axis)];
		addMatvec(block[0], -1.0, in, e_off, out, e_off, term.dofs, ndofs);
		addMatvec(block[2], -1.0, in, pe, out, e_off, term.dofs, ndofs);
		addMatvec(block[1], -1.0, in, h_off, out, h_off, term.dofs, ndofs);
		addMatvec(block[3], -1.0, in, ph, out, h_off, term.dofs, ndofs);
		addMatvec(block[6], 1.0, in, e_off, out, pe, term.dofs, ndofs);
		addMatvec(block[4], -1.0, in, pe, out, pe, term.dofs, ndofs);
		addMatvec(block[7], 1.0, in, h_off, out, ph, term.dofs, ndofs);
		addMatvec(block[5], -1.0, in, ph, out, ph, term.dofs, ndofs);
	}
}

void applyDebyeTerm(const HesthavenDebyeTerm& term, int ndofs, int p_base,
                    const double* in, double* out)
{
	for (int axis = 0; axis < 3; ++axis) {
		const int e_off = axis * ndofs;
		const int p_off = p_base + axis * ndofs;
		for (int dof : term.dofs) {
			const double e = in[e_off + dof];
			const double p = in[p_off + dof];
			out[e_off + dof] += term.e_from_e * e + term.e_from_p * p;
			out[p_off + dof] += term.p_from_e * e + term.p_from_p * p;
		}
	}
}

void applyLorentzTerm(const HesthavenLorentzTerm& term, int ndofs, int p_base, int j_base,
                      const double* in, double* out)
{
	for (int axis = 0; axis < 3; ++axis) {
		const int e_off = axis * ndofs;
		const int p_off = p_base + axis * ndofs;
		const int j_off = j_base + axis * ndofs;
		if (!term.spatial) {
			for (int dof : term.dofs) {
				const double e = in[e_off + dof];
				const double p = in[p_off + dof];
				const double j = in[j_off + dof];
				out[e_off + dof] += term.e_from_j * j;
				out[p_off + dof] += j;
				out[j_off + dof] += term.j_from_e * e + term.j_from_p * p + term.j_from_j * j;
			}
			continue;
		}
		for (int dof : term.dofs) {
			const double j = in[j_off + dof];
			out[e_off + dof] += term.e_from_j * j;
			out[p_off + dof] += j;
		}
		addMatvec(term.j_from_e_block, 1.0, in, e_off, out, j_off, term.dofs, ndofs);
		if (term.has_p_block) {
			addMatvec(term.j_from_p_block, 1.0, in, p_off, out, j_off, term.dofs, ndofs);
		}
		addMatvec(term.j_from_j_block, 1.0, in, j_off, out, j_off, term.dofs, ndofs);
	}
}

} // namespace

void HesthavenAuxiliary::build(mfem::ParFiniteElementSpace& fes, Model& model,
                               const mfem::Array<int>& linear_ids)
{
	ndofs_ = fes.GetNDofs();
	const int dim = fes.GetMesh()->Dimension();
	n_pml_ = computePMLAuxSize(model.getPMLProperties(), ndofs_, dim);
	n_debye_ = model.debyeAuxSize(ndofs_);
	n_lorentz_ = model.lorentzAuxSize(ndofs_);
	pml_needs_a_ = false;
	for (const PMLProperties& props : model.getPMLProperties()) {
		if (props.kappa_max > 1.0 + 1e-12) {
			pml_needs_a_ = true;
		}
	}

	linear_.assign(static_cast<size_t>(linear_ids.Size()), HesthavenLinearMaterial{});
	any_linear_scale_ = false;
	const auto& materials = model.getGeomTagToMaterial();
	for (int slot = 0; slot < linear_ids.Size(); ++slot) {
		const int attr = fes.GetMesh()->GetAttribute(linear_ids[slot]);
		const auto it = materials.find(attr);
		if (it == materials.end()) {
			continue;
		}
		const double eps = it->second.getPermittivity();
		const double mu = it->second.getPermeability();
		const double sigma = it->second.getConductivity();
		HesthavenLinearMaterial mat;
		mat.inv_eps = (eps == 1.0) ? 1.0 : 1.0 / eps;
		mat.inv_mu = (mu == 1.0) ? 1.0 : 1.0 / mu;
		mat.sigma_over_eps = (sigma == 0.0) ? 0.0 : sigma / eps;
		linear_[static_cast<size_t>(slot)] = mat;
		any_linear_scale_ = any_linear_scale_ || mat.active();
	}

	debye_.clear();
	lorentz_.clear();
	pml_.clear();
	if (model.hasPML()) {
		const PMLProfileData* profiles = model.getPMLProfileData();
		if (profiles == nullptr) {
			throw std::runtime_error(
				"Hesthaven SC-PML requires initialized PML profile data.");
		}
	}

	for (int elem = 0; elem < fes.GetNE(); ++elem) {
		const int attr = fes.GetMesh()->GetAttribute(elem);
		if (model.isPMLAttribute(attr)) {
			pml_.push_back(buildPmlElement(
				fes, model, elem, pml_needs_a_, *model.getPMLProfileData()));
		}
		if (const DebyeProperties* pole = model.findDebye(attr)) {
			HesthavenDebyeTerm term;
			term.dofs = elementDofs(fes, elem);
			const double eps_d = pole->eps_s - pole->eps_inf;
			const double tau = pole->tau_solver;
			const double eps_inf = pole->eps_inf;
			term.e_from_e = -eps_d / (eps_inf * tau);
			term.e_from_p = 1.0 / (eps_inf * tau);
			term.p_from_e = eps_d / tau;
			term.p_from_p = -1.0 / tau;
			debye_.push_back(std::move(term));
		}
		if (const LorentzProperties* pole = model.findLorentz(attr)) {
			lorentz_.push_back(buildLorentzTerm(fes, elem, *pole));
		}
	}

	if (!pml_.empty()) {
		buildPmlCsr();
	}
}

void HesthavenAuxiliary::buildPmlCsr()
{
	pml_csr_.reset();
	if (pml_.empty() || n_pml_ == 0) {
		return;
	}
	// Full ODE height so AddMult accepts the evolution state vector.
	const int height = 6 * ndofs_ + n_pml_ + n_debye_ + n_lorentz_;
	auto csr = std::make_unique<mfem::SparseMatrix>(height, height);

	// in → out ADE blocks only. Kappa stretch (delta) is out += Δ·out and is applied after SpMV.
	for (const HesthavenPMLElement& term : pml_) {
		for (int axis = 0; axis < 3; ++axis) {
			const int e_off = axis * ndofs_;
			const int h_off = (3 + axis) * ndofs_;
			const int pe = 6 * ndofs_ + axis * ndofs_;
			const int ph = 6 * ndofs_ + (3 + axis) * ndofs_;
			const auto& block = term.blocks[static_cast<size_t>(axis)];
			accumulateDenseBlock(*csr, block[0], -1.0, term.dofs, e_off, e_off);
			accumulateDenseBlock(*csr, block[2], -1.0, term.dofs, e_off, pe);
			accumulateDenseBlock(*csr, block[1], -1.0, term.dofs, h_off, h_off);
			accumulateDenseBlock(*csr, block[3], -1.0, term.dofs, h_off, ph);
			accumulateDenseBlock(*csr, block[6], 1.0, term.dofs, pe, e_off);
			accumulateDenseBlock(*csr, block[4], -1.0, term.dofs, pe, pe);
			accumulateDenseBlock(*csr, block[7], 1.0, term.dofs, ph, h_off);
			accumulateDenseBlock(*csr, block[5], -1.0, term.dofs, ph, ph);
		}
	}

	csr->Finalize();
	pml_csr_ = std::move(csr);
}

const HesthavenLinearMaterial& HesthavenAuxiliary::linearMaterial(int slot) const
{
	if (slot < 0 || slot >= static_cast<int>(linear_.size())) {
		return kVacuumMaterial;
	}
	return linear_[static_cast<size_t>(slot)];
}

void HesthavenAuxiliary::apply(const mfem::Vector& in, mfem::Vector& out) const
{
	if (!hasAux()) {
		return;
	}
	const double* in_data = in.HostRead();
	double* out_data = out.HostReadWrite();
	const int p_debye = 6 * ndofs_ + n_pml_;
	const int p_lorentz = p_debye + n_debye_;
	const int j_lorentz = p_lorentz + 3 * ndofs_;

	if (pml_csr_ != nullptr) {
		// Match applyPmlElement: out += Δ·out on the Maxwell RHS, then ADE SpMV from in.
		if (pml_needs_a_) {
#ifdef SEMBA_DGTD_ENABLE_OPENMP
			#pragma omp parallel for schedule(static) if (pml_.size() > 1)
#endif
			for (int i = 0; i < static_cast<int>(pml_.size()); ++i) {
				const auto& term = pml_[static_cast<size_t>(i)];
				if (!term.needs_a) {
					continue;
				}
				for (int axis = 0; axis < 3; ++axis) {
					const int e_off = axis * ndofs_;
					const int h_off = (3 + axis) * ndofs_;
					addMatvec(term.delta[static_cast<size_t>(axis)], 1.0,
					          out_data, e_off, out_data, e_off, term.dofs, ndofs_);
					addMatvec(term.delta[static_cast<size_t>(axis)], 1.0,
					          out_data, h_off, out_data, h_off, term.dofs, ndofs_);
				}
			}
		}
		pml_csr_->AddMult(in, out);
	} else {
#ifdef SEMBA_DGTD_ENABLE_OPENMP
		#pragma omp parallel for schedule(static) if (pml_.size() > 1)
#endif
		for (int i = 0; i < static_cast<int>(pml_.size()); ++i) {
			applyPmlElement(pml_[static_cast<size_t>(i)], ndofs_, in_data, out_data);
		}
	}
#ifdef SEMBA_DGTD_ENABLE_OPENMP
	#pragma omp parallel for schedule(static) if (debye_.size() > 1)
#endif
	for (int i = 0; i < static_cast<int>(debye_.size()); ++i) {
		applyDebyeTerm(debye_[static_cast<size_t>(i)], ndofs_, p_debye, in_data, out_data);
	}
#ifdef SEMBA_DGTD_ENABLE_OPENMP
	#pragma omp parallel for schedule(static) if (lorentz_.size() > 1)
#endif
	for (int i = 0; i < static_cast<int>(lorentz_.size()); ++i) {
		applyLorentzTerm(lorentz_[static_cast<size_t>(i)], ndofs_, p_lorentz, j_lorentz,
		                 in_data, out_data);
	}
}

} // namespace maxwell
