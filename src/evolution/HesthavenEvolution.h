#pragma once

#include "components/SubMesher.h"

#include "evolution/MaxwellEvolutionMethods.h"
#include "evolution/HesthavenEvolutionMethods.h"

#include <optional>
#include <array>

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
};
#endif

class HesthavenEvolution : public mfem::TimeDependentOperator
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	HesthavenEvolution(mfem::ParFiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& in, mfem::Vector& out) const;
	mfem::MemoryClass GetMemoryClass() const override;

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
	void initGPUBoundaryData();
#endif
	void exchangeFieldData(const mfem::Vector& in, bool for_gpu_mult) const;
	void computeJumps(const mfem::Vector& in, HesthavenFields& jumps) const;

	void evaluateTFSF(HesthavenFields& jumps) const;
	void evaluateTFSF(double* e_jumps, double* h_jumps, int jumps_size) const;
	void addCurvedElementContributions(mfem::Vector& out) const;
	void storeDirectionalMatrices(FiniteElementSpace& subFES, const DynamicMatrix& refInvMass, HesthavenElement&);
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

	Array<ElementId> linearElements_;
	std::map<ElementId, Array<NodeId>> curvedElements_;

	std::set<DynamicMatrix, MatrixCompareLessThan> matrixStorage_;
	std::vector<HesthavenElement> hestElemLinearStorage_;
	std::vector<HesthavenCurvedElement> hestElemCurvedStorage_;
	DynamicMatrix refLIFT_;
	std::optional<Connectivities> connectivity_;
	std::vector<Source::Position> positions_;

	mutable std::array<mfem::ParGridFunction, 3> eOld_;
	mutable std::array<mfem::ParGridFunction, 3> hOld_;

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
#endif

} // namespace maxwell
