#pragma once

#include "evolution/HesthavenEvolutionMethods.h"

#include <array>
#include <memory>
#include <vector>

namespace mfem {
class Vector;
class ParFiniteElementSpace;
class SparseMatrix;
}

namespace maxwell {

class Model;

/// Piecewise-constant scale for one straight element. Vacuum is the identity.
struct HesthavenLinearMaterial {
	double inv_eps = 1.0;
	double inv_mu = 1.0;
	double sigma_over_eps = 0.0;

	bool active() const
	{
		return inv_eps != 1.0 || inv_mu != 1.0 || sigma_over_eps != 0.0;
	}
};

struct HesthavenDebyeTerm {
	std::vector<int> dofs;
	double e_from_e = 0.0;
	double e_from_p = 0.0;
	double p_from_e = 0.0;
	double p_from_p = 0.0;
};

struct HesthavenLorentzTerm {
	std::vector<int> dofs;
	double e_from_j = 0.0;
	double j_from_e = 0.0;
	double j_from_p = 0.0;
	double j_from_j = 0.0;
	bool spatial = false;
	bool has_p_block = false;
	DynamicMatrix j_from_e_block;
	DynamicMatrix j_from_p_block;
	DynamicMatrix j_from_j_block;
};

/// Element blocks matching `DGOperatorFactory::buildSCPMLOperators`.
/// Kinds: 0 bE, 1 bH, 2 cE, 3 cH, 4 dE, 5 dH, 6 invkE, 7 invkH.
struct HesthavenPMLElement {
	std::vector<int> dofs;
	bool needs_a = false;
	std::array<std::array<DynamicMatrix, 8>, 3> blocks;
	std::array<DynamicMatrix, 3> delta;
};

class HesthavenAuxiliary {
public:
	void build(mfem::ParFiniteElementSpace& fes, Model& model, const mfem::Array<int>& linear_ids);

	const HesthavenLinearMaterial& linearMaterial(int slot) const;
	bool anyLinearScale() const { return any_linear_scale_; }
	bool hasAux() const { return !pml_.empty() || !debye_.empty() || !lorentz_.empty(); }

	int ndofs() const { return ndofs_; }
	int nPml() const { return n_pml_; }
	int nDebye() const { return n_debye_; }
	int nLorentz() const { return n_lorentz_; }
	bool pmlNeedsA() const { return pml_needs_a_; }

	const std::vector<HesthavenLinearMaterial>& linearMaterials() const { return linear_; }
	const std::vector<HesthavenDebyeTerm>& debye() const { return debye_; }
	const std::vector<HesthavenLorentzTerm>& lorentz() const { return lorentz_; }
	const std::vector<HesthavenPMLElement>& pml() const { return pml_; }

	/// PML, then Debye, then Lorentz. No-op when none of those are present.
	void apply(const mfem::Vector& in, mfem::Vector& out) const;

private:
	void buildPmlCsr();

	int ndofs_ = 0;
	int n_pml_ = 0;
	int n_debye_ = 0;
	int n_lorentz_ = 0;
	bool pml_needs_a_ = false;
	bool any_linear_scale_ = false;
	std::vector<HesthavenLinearMaterial> linear_;
	std::vector<HesthavenDebyeTerm> debye_;
	std::vector<HesthavenLorentzTerm> lorentz_;
	std::vector<HesthavenPMLElement> pml_;
	/// CPU Mult: merged element PML blocks (Global-style SpMV). Dense pml_ kept for CUDA.
	std::unique_ptr<mfem::SparseMatrix> pml_csr_;
};

} // namespace maxwell
