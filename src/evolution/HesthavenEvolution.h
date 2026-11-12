#pragma once

#include "components/SubMesher.h"

#include "evolution/MaxwellEvolutionMethods.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "evolution/HesthavenAuxiliary.h"

#include <optional>
#include <array>
#include <memory>
#include <unordered_set>

namespace maxwell {

class SourcesManager;
using FieldGridFuncs = std::array<std::array<mfem::GridFunction, 3>, 2>;

	using straightElemList = mfem::Array<ElementId>;
	using curvedElemMap = std::map<ElementId, mfem::Array<NodeId>>;
	std::pair<straightElemList, curvedElemMap> initCurvedAndLinearElementsLists(const mfem::ParFiniteElementSpace& fes, const std::vector<Source::Position>& curved_pos);

#ifdef SEMBA_DGTD_ENABLE_CUDA
/// Device-side layout for Hesthaven element-local operators.
struct HesthavenGPUData {
	bool initialized = false;
	int ndofs = 0;
	int n_elements = 0;
	int n_matrices = 0;
	int jumps_size = 0;
	int lift_rows = 0;
	int lift_cols = 0;
	int ndof_el = 0;
	int max_flux_size = 0;
	int team_size = 0;
	int workspace_stride = 0;

	mfem::Vector d_matrices;
	mfem::Array<int> d_matrix_offsets;
	mfem::Array<int> d_matrix_rows;
	mfem::Array<int> d_matrix_cols;

	mfem::Vector d_ref_lift;

	mfem::Array<int> d_elem_ids;
	mfem::Array<int> d_elem_dofs;
	mfem::Array<int> d_dir_matrix_id;
	mfem::Array<int> d_flux_size;
	mfem::Vector d_normals;
	mfem::Vector d_fscale;

	mfem::Array<int> d_jump_minus;
	mfem::Array<int> d_jump_plus;
	mfem::Array<uint8_t> d_jump_plus_is_nbr;

	mfem::Vector d_jumps_e;
	mfem::Vector d_jumps_h;

	mfem::Vector d_elem_out_e;
	mfem::Vector d_elem_out_h;
	mfem::Vector d_workspace;

	int n_bc_true = 0;
	int n_bc_int = 0;
	int n_tfsf = 0;

	mfem::Array<int> d_bc_true_jump_out;
	mfem::Array<int> d_bc_true_vmap_in;
	mfem::Array<int> d_bc_true_comp;
	mfem::Vector d_bc_true_e_coeff;
	mfem::Vector d_bc_true_h_coeff;

	mfem::Array<int> d_bc_int_jump_out1;
	mfem::Array<int> d_bc_int_jump_out2;
	mfem::Array<int> d_bc_int_vmap_in1;
	mfem::Array<int> d_bc_int_vmap_in2;
	mfem::Array<int> d_bc_int_comp;
	mfem::Vector d_bc_int_e_coeff;
	mfem::Vector d_bc_int_h_coeff;

	mfem::Array<int> d_tfsf_jump_sf;
	mfem::Array<int> d_tfsf_jump_tf;
	mfem::Array<int> d_tfsf_src_dof;

	bool has_tfsf = false;

	mfem::Vector d_inv_eps;
	mfem::Vector d_inv_mu;
	mfem::Vector d_sigma_over_eps;

	int n_debye_elem = 0;
	int debye_p_base = 0;
	mfem::Array<int> d_debye_dofs;
	mfem::Vector d_debye_coeff;

	int n_lorentz_scalar = 0;
	int n_lorentz_spatial = 0;
	int lorentz_p_base = 0;
	int lorentz_j_base = 0;
	mfem::Array<int> d_lorentz_dofs;
	mfem::Vector d_lorentz_coeff;
	mfem::Array<int> d_lorentz_sp_dofs;
	mfem::Vector d_lorentz_sp_e_from_j;
	mfem::Array<int> d_lorentz_sp_has_p;
	mfem::Vector d_lorentz_sp_j_from_e;
	mfem::Vector d_lorentz_sp_j_from_p;
	mfem::Vector d_lorentz_sp_j_from_j;

	int n_pml_elem = 0;
	int pml_needs_a = 0;
	mfem::Array<int> d_pml_dofs;
	mfem::Vector d_pml_blocks;
	mfem::Vector d_pml_delta;

	bool curved_initialized = false;
	int n_curved_rows = 0;
	int curved_full_size = 0;
	mfem::Vector d_ext_in;
	mfem::Vector d_curved_y;
	mfem::Array<int> d_curved_row_to_out;
};
#endif

class GlobalEvolution;
class Probes;

class HesthavenEvolution : public mfem::TimeDependentOperator
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	HesthavenEvolution(mfem::ParFiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&,
	                   const Probes* probes = nullptr, double final_time = 0.0);
	~HesthavenEvolution();
	void SetTime(double t) override;
	virtual void Mult(const mfem::Vector& in, mfem::Vector& out) const;
	void ImplicitSolve(double dt, const mfem::Vector& x, mfem::Vector& k) override;
	mfem::MemoryClass GetMemoryClass() const override;

	GlobalEvolution* sgbcEngine() { return sgbc_engine_.get(); }
	const GlobalEvolution* sgbcEngine() const { return sgbc_engine_.get(); }

	const HesthavenElement& getHesthavenElement(const ElementId& e) const { return hestElemLinearStorage_[e]; }

#ifdef SEMBA_DGTD_ENABLE_CUDA
	/// CPU Mult path for GPU-vs-CPU benchmark tests only.
	void benchmarkMultCPU(const mfem::Vector& in, mfem::Vector& out) const { MultCPU(in, out); }
#endif

private:
	void MultCPU(const mfem::Vector& in, mfem::Vector& out) const;
#ifdef SEMBA_DGTD_ENABLE_CUDA
	void MultGPU(const mfem::Vector& in, mfem::Vector& out) const;
	void initGPUData();
	void initGPUCurvedData();
	void initGPUBoundaryData();
#endif
	void exchangeFieldData(const mfem::Vector& in, bool for_gpu_mult) const;
	void promoteAssembledElements();
	void applySGBC(const mfem::Vector& in, mfem::Vector& out, bool inject) const;
	void computeJumps(const mfem::Vector& in, HesthavenFields& jumps) const;
	void buildMergedCurvedMatrix();

	void evaluateTFSF(HesthavenFields& jumps) const;
	void evaluateTFSF(double* e_jumps, double* h_jumps, int jumps_size) const;
	void addCurvedElementContributions(mfem::Vector& out) const;
	void storeDirectionalMatrices(const DynamicMatrix& refInvMass, HesthavenElement&);
	void storeFaceInformation(HesthavenElement&);
	void buildLinearDirGroups();
	void checkForTFSFInCurvedElements();
	void applyBoundaryConditionsToNodes(const BoundaryMaps&, const FieldsInputMaps& in, HesthavenFields& out) const;
	void applyBoundaryConditionsToNodes(const BoundaryMaps&, const FieldsInputMaps& in,
	                                    double* e_jumps, double* h_jumps, int jumps_size) const;
	bool isDoFinCurvedElement(const NodeId& d) const;

	const Eigen::VectorXd applyLIFT(const Eigen::VectorXd& fscale, Eigen::VectorXd& flux) const;

	ParFiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;
	HesthavenAuxiliary auxiliary_;

	Array<ElementId> linearElements_;
	std::map<ElementId, Array<NodeId>> curvedElements_;

	std::set<DynamicMatrix, MatrixCompareLessThan> matrixStorage_;
	std::vector<HesthavenElement> hestElemLinearStorage_;
	/// Slots in hestElemLinearStorage_ that share the same dir[0..2] pointers.
	std::vector<std::vector<int>> linear_dir_groups_;
	std::vector<HesthavenCurvedElement> hestElemCurvedStorage_;
	DynamicMatrix refLIFT_;
	std::optional<Connectivities> connectivity_;
	std::vector<Source::Position> positions_;

	/// Preallocated per-thread scratch for the element-local Mult (no heap
	/// allocation inside MultCPU). Sized to the linear elements' DOF / flux sizes.
	struct MultScratch {
		Eigen::VectorXd eIn[3], hIn[3];
		Eigen::VectorXd ndotdH, ndotdE;
		Eigen::VectorXd fluxH, fluxE, tmpF, tmpF2;
		Eigen::VectorXd liftedH, liftedE, volH, volE, tmpV;
		void ensure(int nd, int nf)
		{
			if (volH.size() != nd) {
				for (int d = X; d <= Z; ++d) {
					eIn[d].resize(nd);
					hIn[d].resize(nd);
				}
				liftedH.resize(nd);
				liftedE.resize(nd);
				volH.resize(nd);
				volE.resize(nd);
				tmpV.resize(nd);
			}
			if (fluxH.size() != nf) {
				ndotdH.resize(nf);
				ndotdE.resize(nf);
				fluxH.resize(nf);
				fluxE.resize(nf);
				tmpF.resize(nf);
				tmpF2.resize(nf);
			}
		}
	};
	mutable std::vector<MultScratch> mult_scratch_;
	/// Preallocated jump buffer (element-major, connectivity_->global.size()).
	mutable std::unique_ptr<HesthavenFields> mult_jumps_;

	mutable std::array<mfem::ParGridFunction, 3> eOld_;
	mutable std::array<mfem::ParGridFunction, 3> hOld_;

	std::unique_ptr<GlobalEvolution> sgbc_engine_;
	/// Owns options for nested SGBC GlobalEvolution (must outlive sgbc_engine_).
	std::unique_ptr<EvolutionOptions> sgbc_host_opts_;

	std::unique_ptr<mfem::SparseMatrix> curved_merged_matrix_;
	mfem::Array<int> curved_merged_row_to_out_;

	mutable mfem::Vector implicit_rhs_;
	mutable mfem::Vector implicit_tmp_;
	mutable double implicit_cached_dt_{-1.0};

#ifdef SEMBA_DGTD_ENABLE_CUDA
	mutable HesthavenGPUData gpu_;
	bool has_tfsf_gpu_ = false;
#endif
};

#ifdef SEMBA_DGTD_ENABLE_CUDA
void hesthaven_compute_jumps_gpu(HesthavenGPUData& gpu,
                                 const mfem::Vector& in,
                                 const std::array<mfem::ParGridFunction, 3>& eOld,
                                 const std::array<mfem::ParGridFunction, 3>& hOld,
                                 int ndofs);

void hesthaven_mult_gpu(HesthavenGPUData& gpu,
                        double alpha,
                        const mfem::Vector& in,
                        const std::array<mfem::ParGridFunction, 3>& eOld,
                        const std::array<mfem::ParGridFunction, 3>& hOld,
                        mfem::Vector& out,
                        int ndofs);

void hesthaven_apply_bc_gpu(HesthavenGPUData& gpu,
                            const mfem::Vector& in,
                            int ndofs);

void hesthaven_scatter_tfsf_to_jumps_gpu(HesthavenGPUData& gpu,
                                         const FieldGridFuncs& tfsf_fields,
                                         int jumps_size);

void hesthaven_apply_aux_gpu(HesthavenGPUData& gpu,
                             const mfem::Vector& in,
                             mfem::Vector& out,
                             int ndofs);

void hesthaven_add_curved_gpu(HesthavenGPUData& gpu,
                              mfem::SparseMatrix& curved_merged,
                              const std::array<mfem::ParGridFunction, 3>& eOld,
                              const std::array<mfem::ParGridFunction, 3>& hOld,
                              mfem::Vector& out,
                              int ndofs,
                              int nbr_size);
#endif

} // namespace maxwell
