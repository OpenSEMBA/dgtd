#include "HesthavenEvolution.h"

#include "components/DGOperatorFactory.h"
#include "solver/SourcesManager.h"
#include "GlobalEvolution.h"

#include "components/SCPMLLayout.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#ifdef SEMBA_DGTD_ENABLE_OPENMP
#include <omp.h>
#endif

namespace maxwell {

#ifdef SEMBA_DGTD_ENABLE_CUDA
namespace {

bool hesthavenGpuMultViable()
{
	return mfem::Device::Allows(mfem::Backend::CUDA);
}

void hesthaven_sync_gpu_static_data(HesthavenGPUData& gpu)
{
	// HostWrite() fills above; push to device only (ReadWrite would pull stale device data).
	auto sync_vec = [](mfem::Vector& v) {
		if (v.Size() > 0) { (void)v.Read(); }
	};
	auto sync_arr = [](auto& a) {
		if (a.Size() > 0) { (void)a.Read(); }
	};

	sync_vec(gpu.d_matrices);
	sync_vec(gpu.d_ref_lift);
	sync_arr(gpu.d_matrix_offsets);
	sync_arr(gpu.d_matrix_rows);
	sync_arr(gpu.d_matrix_cols);
	sync_arr(gpu.d_elem_ids);
	sync_arr(gpu.d_elem_dofs);
	sync_arr(gpu.d_dir_matrix_id);
	sync_arr(gpu.d_flux_size);
	sync_vec(gpu.d_normals);
	sync_vec(gpu.d_fscale);
	sync_arr(gpu.d_jump_minus);
	sync_arr(gpu.d_jump_plus);
	sync_arr(gpu.d_jump_plus_is_nbr);
	sync_arr(gpu.d_bc_true_jump_out);
	sync_arr(gpu.d_bc_true_vmap_in);
	sync_arr(gpu.d_bc_true_comp);
	sync_vec(gpu.d_bc_true_e_coeff);
	sync_vec(gpu.d_bc_true_h_coeff);
	sync_arr(gpu.d_bc_int_jump_out1);
	sync_arr(gpu.d_bc_int_jump_out2);
	sync_arr(gpu.d_bc_int_vmap_in1);
	sync_arr(gpu.d_bc_int_vmap_in2);
	sync_arr(gpu.d_bc_int_comp);
	sync_vec(gpu.d_bc_int_e_coeff);
	sync_vec(gpu.d_bc_int_h_coeff);
	sync_arr(gpu.d_tfsf_jump_sf);
	sync_arr(gpu.d_tfsf_jump_tf);
	sync_arr(gpu.d_tfsf_src_dof);
	sync_vec(gpu.d_inv_eps);
	sync_vec(gpu.d_inv_mu);
	sync_vec(gpu.d_sigma_over_eps);
	sync_arr(gpu.d_debye_dofs);
	sync_vec(gpu.d_debye_coeff);
	sync_arr(gpu.d_lorentz_dofs);
	sync_vec(gpu.d_lorentz_coeff);
	sync_arr(gpu.d_lorentz_sp_dofs);
	sync_vec(gpu.d_lorentz_sp_e_from_j);
	sync_arr(gpu.d_lorentz_sp_has_p);
	sync_vec(gpu.d_lorentz_sp_j_from_e);
	sync_vec(gpu.d_lorentz_sp_j_from_p);
	sync_vec(gpu.d_lorentz_sp_j_from_j);
	sync_arr(gpu.d_pml_dofs);
	sync_vec(gpu.d_pml_blocks);
	sync_vec(gpu.d_pml_delta);
}

void pack_element_dofs(const std::vector<int>& dofs, int ndof_el, int* dst)
{
	if (static_cast<int>(dofs.size()) != ndof_el) {
		throw std::runtime_error(
			"Hesthaven GPU auxiliary path requires a uniform element dof count.");
	}
	for (int i = 0; i < ndof_el; ++i) {
		dst[i] = dofs[static_cast<size_t>(i)];
	}
}

void packHesthavenAuxGpu(const HesthavenAuxiliary& aux, HesthavenGPUData& gpu)
{
	const int ndof = gpu.ndof_el;
	const int n2 = ndof * ndof;
	if (aux.anyLinearScale()) {
		const int n = static_cast<int>(aux.linearMaterials().size());
		gpu.d_inv_eps.SetSize(n);
		gpu.d_inv_mu.SetSize(n);
		gpu.d_sigma_over_eps.SetSize(n);
		gpu.d_inv_eps.UseDevice(true);
		gpu.d_inv_mu.UseDevice(true);
		gpu.d_sigma_over_eps.UseDevice(true);
		double* inv_eps = gpu.d_inv_eps.HostWrite();
		double* inv_mu = gpu.d_inv_mu.HostWrite();
		double* sigma = gpu.d_sigma_over_eps.HostWrite();
		for (int e = 0; e < n; ++e) {
			const HesthavenLinearMaterial& mat = aux.linearMaterials()[static_cast<size_t>(e)];
			inv_eps[e] = mat.inv_eps;
			inv_mu[e] = mat.inv_mu;
			sigma[e] = mat.sigma_over_eps;
		}
	}

	gpu.n_debye_elem = static_cast<int>(aux.debye().size());
	gpu.debye_p_base = 6 * aux.ndofs() + aux.nPml();
	if (gpu.n_debye_elem > 0) {
		gpu.d_debye_dofs.SetSize(gpu.n_debye_elem * ndof);
		gpu.d_debye_dofs.GetMemory().UseDevice(true);
		gpu.d_debye_coeff.SetSize(gpu.n_debye_elem * 4);
		gpu.d_debye_coeff.UseDevice(true);
		int* dofs = gpu.d_debye_dofs.HostWrite();
		double* coeff = gpu.d_debye_coeff.HostWrite();
		for (int e = 0; e < gpu.n_debye_elem; ++e) {
			const HesthavenDebyeTerm& term = aux.debye()[static_cast<size_t>(e)];
			pack_element_dofs(term.dofs, ndof, dofs + e * ndof);
			coeff[4 * e + 0] = term.e_from_e;
			coeff[4 * e + 1] = term.e_from_p;
			coeff[4 * e + 2] = term.p_from_e;
			coeff[4 * e + 3] = term.p_from_p;
		}
	}

	gpu.lorentz_p_base = gpu.debye_p_base + aux.nDebye();
	gpu.lorentz_j_base = gpu.lorentz_p_base + 3 * aux.ndofs();
	int n_scalar = 0;
	int n_spatial = 0;
	for (const HesthavenLorentzTerm& term : aux.lorentz()) {
		if (term.spatial) {
			++n_spatial;
		} else {
			++n_scalar;
		}
	}
	gpu.n_lorentz_scalar = n_scalar;
	gpu.n_lorentz_spatial = n_spatial;
	if (n_scalar > 0) {
		gpu.d_lorentz_dofs.SetSize(n_scalar * ndof);
		gpu.d_lorentz_dofs.GetMemory().UseDevice(true);
		gpu.d_lorentz_coeff.SetSize(n_scalar * 4);
		gpu.d_lorentz_coeff.UseDevice(true);
		int* dofs = gpu.d_lorentz_dofs.HostWrite();
		double* coeff = gpu.d_lorentz_coeff.HostWrite();
		int slot = 0;
		for (const HesthavenLorentzTerm& term : aux.lorentz()) {
			if (term.spatial) {
				continue;
			}
			pack_element_dofs(term.dofs, ndof, dofs + slot * ndof);
			coeff[4 * slot + 0] = term.e_from_j;
			coeff[4 * slot + 1] = term.j_from_e;
			coeff[4 * slot + 2] = term.j_from_p;
			coeff[4 * slot + 3] = term.j_from_j;
			++slot;
		}
	}
	if (n_spatial > 0) {
		gpu.d_lorentz_sp_dofs.SetSize(n_spatial * ndof);
		gpu.d_lorentz_sp_dofs.GetMemory().UseDevice(true);
		gpu.d_lorentz_sp_e_from_j.SetSize(n_spatial);
		gpu.d_lorentz_sp_e_from_j.UseDevice(true);
		gpu.d_lorentz_sp_has_p.SetSize(n_spatial);
		gpu.d_lorentz_sp_has_p.GetMemory().UseDevice(true);
		gpu.d_lorentz_sp_j_from_e.SetSize(n_spatial * n2);
		gpu.d_lorentz_sp_j_from_p.SetSize(n_spatial * n2);
		gpu.d_lorentz_sp_j_from_j.SetSize(n_spatial * n2);
		gpu.d_lorentz_sp_j_from_e.UseDevice(true);
		gpu.d_lorentz_sp_j_from_p.UseDevice(true);
		gpu.d_lorentz_sp_j_from_j.UseDevice(true);
		int* dofs = gpu.d_lorentz_sp_dofs.HostWrite();
		double* e_from_j = gpu.d_lorentz_sp_e_from_j.HostWrite();
		int* has_p = gpu.d_lorentz_sp_has_p.HostWrite();
		double* je = gpu.d_lorentz_sp_j_from_e.HostWrite();
		double* jp = gpu.d_lorentz_sp_j_from_p.HostWrite();
		double* jj = gpu.d_lorentz_sp_j_from_j.HostWrite();
		int slot = 0;
		for (const HesthavenLorentzTerm& term : aux.lorentz()) {
			if (!term.spatial) {
				continue;
			}
			pack_element_dofs(term.dofs, ndof, dofs + slot * ndof);
			e_from_j[slot] = term.e_from_j;
			has_p[slot] = term.has_p_block ? 1 : 0;
			std::memcpy(je + slot * n2, term.j_from_e_block.data(),
			            static_cast<size_t>(n2) * sizeof(double));
			if (term.has_p_block) {
				std::memcpy(jp + slot * n2, term.j_from_p_block.data(),
				            static_cast<size_t>(n2) * sizeof(double));
			} else {
				std::fill(jp + slot * n2, jp + (slot + 1) * n2, 0.0);
			}
			std::memcpy(jj + slot * n2, term.j_from_j_block.data(),
			            static_cast<size_t>(n2) * sizeof(double));
			++slot;
		}
	}

	gpu.n_pml_elem = static_cast<int>(aux.pml().size());
	gpu.pml_needs_a = aux.pmlNeedsA() ? 1 : 0;
	if (gpu.n_pml_elem > 0) {
		gpu.d_pml_dofs.SetSize(gpu.n_pml_elem * ndof);
		gpu.d_pml_dofs.GetMemory().UseDevice(true);
		gpu.d_pml_blocks.SetSize(gpu.n_pml_elem * 3 * 8 * n2);
		gpu.d_pml_blocks.UseDevice(true);
		int* dofs = gpu.d_pml_dofs.HostWrite();
		double* blocks = gpu.d_pml_blocks.HostWrite();
		double* delta = nullptr;
		if (gpu.pml_needs_a) {
			gpu.d_pml_delta.SetSize(gpu.n_pml_elem * 3 * n2);
			gpu.d_pml_delta.UseDevice(true);
			delta = gpu.d_pml_delta.HostWrite();
		}
		for (int e = 0; e < gpu.n_pml_elem; ++e) {
			const HesthavenPMLElement& term = aux.pml()[static_cast<size_t>(e)];
			pack_element_dofs(term.dofs, ndof, dofs + e * ndof);
			for (int axis = 0; axis < 3; ++axis) {
				for (int kind = 0; kind < 8; ++kind) {
					const DynamicMatrix& block = term.blocks[static_cast<size_t>(axis)][static_cast<size_t>(kind)];
					const int offset = ((e * 3 + axis) * 8 + kind) * n2;
					std::memcpy(blocks + offset, block.data(),
					            static_cast<size_t>(n2) * sizeof(double));
				}
				if (delta != nullptr) {
					const int offset = (e * 3 + axis) * n2;
					std::memcpy(delta + offset, term.delta[static_cast<size_t>(axis)].data(),
					            static_cast<size_t>(n2) * sizeof(double));
				}
			}
		}
	}
}

} // namespace
#endif

using MatricesSet = std::set<DynamicMatrix, MatrixCompareLessThan>;

int hesthavenStateSize(const ParFiniteElementSpace& fes, const Model& model)
{
	const int ndofs = fes.GetNDofs();
	return 6 * ndofs
		+ computePMLAuxSize(model.getPMLProperties(), ndofs, fes.GetMesh()->Dimension())
		+ model.debyeAuxSize(ndofs)
		+ model.lorentzAuxSize(ndofs);
}

static DynamicMatrix assembleInverseMassMatrix(FiniteElementSpace& fes)
{
	BilinearForm bf(&fes);
	ConstantCoefficient one(1.0);
	bf.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	bf.Assemble();
	bf.Finalize();
	
	auto dense = bf.SpMat().ToDenseMatrix();
	const auto res = toEigen(*dense);
	delete dense;
	return res;
}

Mesh getRefMeshForGeomType(const Element::Type elType, const int dimension)
{
	switch (dimension) {
	case 2:
		switch (elType) {
		case Element::Type::TRIANGLE:
			return Mesh::MakeCartesian2D(1, 1, elType, true, 2.0, 2.0);
		case Element::Type::QUADRILATERAL:
			return Mesh::MakeCartesian2D(1, 1, elType);
		default:
			throw std::runtime_error("Incorrect Element Type for dimension 2 mesh.");
		}
	case 3:
		return buildHesthavenRefTetrahedra();
	default:
		throw std::runtime_error("Hesthaven Evolution Operator only supports dimensions 2 or 3.");
	}
}

DynamicMatrix assembleHesthavenRefElemInvMassMatrix(const Element::Type elType, const int order, const int dimension)
{
	auto m{ getRefMeshForGeomType(elType, dimension) };
	auto fec{ L2_FECollection(order, dimension, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec) };
	auto mass_mat{ assembleInverseMassMatrix(fes) };

	return getElemMassMatrixFromGlobal(0, mass_mat, elType);
}

DynamicMatrix assembleHesthavenRefElemEmat(const Element::Type elType, const int order, const int dimension)
{
	auto m{ getRefMeshForGeomType(elType, dimension) };
	m.SetAttribute(0, hesthavenMeshingTag);
	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);
	auto sm{ SubMesh::CreateFromDomain(m, elementMarker) };
	auto fec{ L2_FECollection(order, dimension, BasisType::GaussLobatto) };
	FiniteElementSpace subFES(&sm, &fec);

	auto boundary_markers = assembleBoundaryMarkers(subFES);

	sm.bdr_attributes.SetSize(boundary_markers.size());
	for (auto f= 0; f < subFES.GetNF(); f++) {
		sm.bdr_attributes[f] = f + 1;
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	DynamicMatrix emat = assembleEmat(subFES, boundary_markers);
	if (dimension == 3 && elType == Element::Type::TETRAHEDRON)
	{
		emat *= 2.0;
	}

	return emat;
}

void initNormalVectors(HesthavenElement& hestElem, const int size)
{
	for (auto d{ X }; d <= Z; d++) {
		hestElem.normals[d].resize(size);
		hestElem.normals[d].setZero();
	}
}

void initFscale(HesthavenElement& hestElem, const int size)
{
	hestElem.fscale.resize(size);
	hestElem.fscale.setZero();
}

void HesthavenEvolution::evaluateTFSF(HesthavenFields& out) const
{
	std::array<std::array<double, 3>, 2> fields;
	const auto& mapBSF = connectivity_->boundary.TFSF.mapBSF;
	const auto& mapBTF = connectivity_->boundary.TFSF.mapBTF;
	const auto& vmapBSF = connectivity_->boundary.TFSF.vmapBSF;
	const auto& vmapBTF = connectivity_->boundary.TFSF.vmapBTF;
	for (const auto& source : srcmngr_.sources) {
		auto tf = dynamic_cast<TotalField*>(source.get());
		if (tf == nullptr) {
			continue;
		}

		for (int m = 0; m < mapBSF.size(); m++) {
			for (int v = 0; v < mapBSF[m].size(); v++) {
				for (int d : { X, Y, Z }) {
					fields[E][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), E, d);
					fields[H][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), H, d);
					// Scattered side subtracts the incident field; total side adds it.
					out.e_[d][mapBSF[m][v]] -= fields[E][d];
					out.h_[d][mapBSF[m][v]] -= fields[H][d];
					out.e_[d][mapBTF[m][v]] += fields[E][d];
					out.h_[d][mapBTF[m][v]] += fields[H][d];
				}
			}
		}
	}
}

void HesthavenEvolution::evaluateTFSF(double* e_jumps, double* h_jumps, int jumps_size) const
{
	std::array<std::array<double, 3>, 2> fields;
	const auto& mapBSF = connectivity_->boundary.TFSF.mapBSF;
	const auto& mapBTF = connectivity_->boundary.TFSF.mapBTF;
	const auto& vmapBSF = connectivity_->boundary.TFSF.vmapBSF;
	const auto& vmapBTF = connectivity_->boundary.TFSF.vmapBTF;
	for (const auto& source : srcmngr_.sources) {
		auto tf = dynamic_cast<TotalField*>(source.get());
		if (tf == nullptr) {
			continue;
		}

		for (int m = 0; m < mapBSF.size(); m++) {
			for (int v = 0; v < mapBSF[m].size(); v++) {
				for (int d : { X, Y, Z }) {
					fields[E][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), E, d);
					fields[H][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), H, d);
					e_jumps[d * jumps_size + mapBSF[m][v]] -= fields[E][d];
					h_jumps[d * jumps_size + mapBSF[m][v]] -= fields[H][d];
					e_jumps[d * jumps_size + mapBTF[m][v]] += fields[E][d];
					h_jumps[d * jumps_size + mapBTF[m][v]] += fields[H][d];
				}
			}
		}
	}

}

void HesthavenEvolution::buildMergedCurvedMatrix()
{
	curved_merged_matrix_.reset();
	curved_merged_row_to_out_.SetSize(0);
	if (hestElemCurvedStorage_.empty()) {
		return;
	}

	const int ndofs = fes_.GetNDofs();
	const int nbrDofs = fes_.num_face_nbr_dofs;
	const int blockSize = ndofs + nbrDofs;
	const int ncomp = numberOfFieldComponents * numberOfMaxDimensions;
	const int fullSize = ncomp * blockSize;

	int n_rows = 0;
	for (const auto& curved : hestElemCurvedStorage_) {
		n_rows += curved.dofs.Size() * ncomp;
	}
	if (n_rows == 0) {
		return;
	}

	curved_merged_matrix_ = std::make_unique<SparseMatrix>(n_rows, fullSize);
	curved_merged_row_to_out_.SetSize(n_rows);

	int row = 0;
	for (const auto& curved : hestElemCurvedStorage_) {
		for (int ii = 0; ii < curved.dofs.Size(); ++ii) {
			const int dof = curved.dofs[ii];
			for (int c = 0; c < ncomp; ++c) {
				const int global_row = c * ndofs + dof;
				Array<int> cols;
				Vector vals;
				curved.matrix.GetRow(global_row, cols, vals);
				curved_merged_matrix_->SetRow(row, cols, vals);
				curved_merged_row_to_out_[row] = global_row;
				++row;
			}
		}
	}
	curved_merged_matrix_->Finalize();
}

void HesthavenEvolution::addCurvedElementContributions(Vector& out) const
{
	if (!curved_merged_matrix_ || curved_merged_matrix_->Height() == 0) {
		return;
	}

	const int ndofs = fes_.GetNDofs();
	const int nbrDofs = fes_.num_face_nbr_dofs;

#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (gpu_.curved_initialized &&
	    mfem::Device::Allows(mfem::Backend::CUDA) && out.UseDevice()) {
		hesthaven_add_curved_gpu(gpu_, *curved_merged_matrix_, eOld_, hOld_, out,
		                         ndofs, nbrDofs);
		return;
	}
#endif

	const int blockSize = ndofs + nbrDofs;
	const int fullSize = numberOfFieldComponents * numberOfMaxDimensions * blockSize;
	const int n_rows = curved_merged_matrix_->Height();

	Vector ext_in(fullSize);
	ext_in.UseDevice(false);
	ext_in = 0.0;
	for (int d = X; d <= Z; ++d) {
		ext_in.SetVector(eOld_[d],       d      * blockSize);
		ext_in.SetVector(hOld_[d],  (3 + d) * blockSize);
		if (nbrDofs > 0) {
			ext_in.SetVector(eOld_[d].FaceNbrData(),      d      * blockSize + ndofs);
			ext_in.SetVector(hOld_[d].FaceNbrData(), (3 + d) * blockSize + ndofs);
		}
	}

	Vector curved_y(n_rows);
	curved_y.UseDevice(false);
	curved_y = 0.0;
	curved_merged_matrix_->AddMult(ext_in, curved_y);

#ifdef SEMBA_DGTD_ENABLE_CUDA
	const bool out_is_device =
	    mfem::Device::Allows(mfem::Backend::CUDA) && out.UseDevice();
	if (out_is_device) {
		(void)out.HostRead();
	}
#endif

	double* out_data = out.HostReadWrite();
	const double* y = curved_y.HostRead();
	for (int row = 0; row < n_rows; ++row) {
		out_data[curved_merged_row_to_out_[row]] += y[row];
	}

#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (out_is_device) {
		(void)out.Read();
	}
#endif
}

const Eigen::VectorXd HesthavenEvolution::applyLIFT(const Eigen::VectorXd& fscale, Eigen::VectorXd& flux) const
{
	for (int i = 0; i < flux.size(); i++) {
		flux[i] *= fscale[i] / 2.0;
	}
	return this->refLIFT_ * flux;
}

double getReferenceVolume(const Element::Type geom)
{
	switch (geom) {
	case Element::Type::TRIANGLE:
		return 2.0; //Hesthaven definition (-1,-1), (-1,1), (1,-1)
	case Element::Type::QUADRILATERAL:
		return 4.0; //Assuming x,y (-1, 1)
	case Element::Type::TETRAHEDRON:
		return 8.0 / 6.0; //Hesthaven definition (-1,-1,-1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)
	case Element::Type::HEXAHEDRON:
		return 8.0; //Assuming x,y,z (-1, 1)
	default:
		throw std::runtime_error("Unsupported geometry for reference volume.");
	}
}

void HesthavenEvolution::storeDirectionalMatrices(const DynamicMatrix& refInvMass, HesthavenElement& hestElem)
{
	// Affine / linear elements: assemble DerivativeIntegrator on the parent element
	// (same integrator as DGOperatorFactory::buildDerivativeSubOperator) — no SubMesh.
	const int e = hestElem.id;
	const FiniteElement& fe = *fes_.GetFE(e);
	ElementTransformation& T = *fes_.GetElementTransformation(e);
	const int dim = fes_.GetMesh()->Dimension();
	const int nd = fe.GetDof();
	ConstantCoefficient one(1.0);
	const double scale = getReferenceVolume(hestElem.type) / hestElem.vol;

	for (int d = X; d <= Z; d++) {
		DenseMatrix elmat(nd);
		elmat = 0.0;
		if (d < dim) {
			DerivativeIntegrator integ(one, d);
			integ.AssembleElementMatrix(fe, T, elmat);
		}
		DynamicMatrix dirMat = refInvMass * toEigen(elmat) * scale;
		StorageIterator it = matrixStorage_.find(dirMat);
		if (it == matrixStorage_.end()) {
			matrixStorage_.insert(dirMat);
			it = matrixStorage_.find(dirMat);
			hestElem.dir[d] = &(*it);
		} else {
			hestElem.dir[d] = &(*it);
		}
	}
}

void HesthavenEvolution::buildLinearDirGroups()
{
	linear_dir_groups_.clear();
	using Key = std::tuple<const DynamicMatrix*, const DynamicMatrix*, const DynamicMatrix*>;
	std::map<Key, std::vector<int>> groups;
	for (int e = 0; e < static_cast<int>(hestElemLinearStorage_.size()); ++e) {
		const auto& he = hestElemLinearStorage_[static_cast<size_t>(e)];
		groups[Key{he.dir[0], he.dir[1], he.dir[2]}].push_back(e);
	}
	linear_dir_groups_.reserve(groups.size());
	for (auto& kv : groups) {
		linear_dir_groups_.push_back(std::move(kv.second));
	}
	if (Mpi::WorldRank() == 0 && !linear_dir_groups_.empty()) {
		std::size_t max_g = 0;
		for (const auto& g : linear_dir_groups_) {
			max_g = std::max(max_g, g.size());
		}
		std::cout << "Hesthaven Mult: " << linear_dir_groups_.size()
		          << " dir-groups (max size " << max_g << ") for batched volume GEMM\n"
		          << std::flush;
	}
}

void HesthavenEvolution::storeFaceInformation(HesthavenElement& hestElem)
{
	Mesh& mesh = *fes_.GetMesh();
	const int e = hestElem.id;
	const int dim = mesh.Dimension();

	Array<int> faces, orients;
	if (dim == 2) {
		mesh.GetElementEdges(e, faces, orients);
	} else {
		mesh.GetElementFaces(e, faces, orients);
	}
	const int numFaces = faces.Size();
	const int numNodesAtFace = (dim == 2)
		? fes_.FEColl()->GetOrder() + 1
		: getFaceNodeNumByGeomType(fes_);

	initNormalVectors(hestElem, numFaces * numNodesAtFace);
	initFscale(hestElem, numFaces * numNodesAtFace);

	const double J = fes_.GetElementTransformation(e)->Weight();

	for (int f = 0; f < numFaces; f++) {
		const int gf = faces[f];
		ElementTransformation* faceTrans = (dim == 2)
			? mesh.GetEdgeTransformation(gf)
			: mesh.GetFaceTransformation(gf);
		const auto* ir = &IntRules.Get(
			faceTrans->GetGeometryType(),
			faceTrans->OrderW() + 2 * faceTrans->Order());
		faceTrans->SetIntPoint(&ir->IntPoint(0));

		Vector normal(dim);
		CalcOrtho(faceTrans->Jacobian(), normal);

		// Parent-mesh edge/face orientation is global; flip so n is outward for element e
		// (MFEM FaceElementTransformations: normal points Elem1 → Elem2).
		FaceElementTransformations* fet = mesh.GetFaceElementTransformations(gf);
		if (fet != nullptr && fet->Elem1No != e) {
			normal *= -1.0;
		}

		const double sJ = faceTrans->Weight();
		for (int b = 0; b < numNodesAtFace; b++) {
			hestElem.normals[X][f * numNodesAtFace + b] = normal[0] / sJ;
			hestElem.fscale[f * numNodesAtFace + b] = sJ * 2.0 / J;
			if (dim >= 2) {
				hestElem.normals[Y][f * numNodesAtFace + b] = normal[1] / sJ;
			}
			if (dim == 3) {
				hestElem.normals[Z][f * numNodesAtFace + b] = normal[2] / sJ;
			}
		}
	}
}

std::pair<Array<ElementId>,std::map<ElementId,Array<NodeId>>> initCurvedAndLinearElementsLists(const ParFiniteElementSpace& fes, const std::vector<Source::Position>& curved_pos)
{
	Mesh mesh_p1(*fes.GetMesh());
	FiniteElementSpace fes_p1(&mesh_p1, fes.FEColl());

	mesh_p1.SetCurvature(1);

	const auto& pos_cur = curved_pos;
	const auto pos_lin = buildDoFPositions(fes_p1);

	double tol{ 1e-5 };
	std::pair<Array<ElementId>, std::map<ElementId, Array<NodeId>>> res;
	for (int e = 0; e < mesh_p1.GetNE(); e++) {
		Array<int> elemdofs_p1, elemdofs_p2;
		fes_p1.GetElementDofs(e, elemdofs_p1);
		fes   .GetElementDofs(e, elemdofs_p2);
		MFEM_ASSERT(elemdofs_p1, elemdofs_p2);
		auto isCurved = false;
		for (auto d= 0; d < elemdofs_p1.Size(); d++) {
			if (std::abs(pos_lin[elemdofs_p1[d]][0] - pos_cur[elemdofs_p2[d]][0]) > tol ||
				std::abs(pos_lin[elemdofs_p1[d]][1] - pos_cur[elemdofs_p2[d]][1]) > tol ||
				std::abs(pos_lin[elemdofs_p1[d]][2] - pos_cur[elemdofs_p2[d]][2]) > tol) 
			{
				isCurved = true;
			}
		}
		if (isCurved) {
			res.second[e] = elemdofs_p2;
		}
		else {
			res.first.Append(e);
		}
	}
	return res;
}

void HesthavenEvolution::checkForTFSFInCurvedElements()
{
	if (model_.getTotalFieldScatteredFieldToMarker().size()) {
		for (const auto& [k, marker] : model_.getTotalFieldScatteredFieldToMarker()) {
			for (auto b= 0; b < fes_.GetNBE(); b++) {
				if (marker[model_.getMesh().GetBdrAttribute(b) - 1] == 1) {
					auto be_trans{ getFaceElementTransformation(model_.getMesh(), b) };
					if (curvedElements_.find(be_trans->Elem1No) != curvedElements_.end()) {
						throw std::runtime_error("TFSF defined on curved elements is not supported.");
					}
					if (be_trans->Elem2No != -1) {
						if (curvedElements_.find(be_trans->Elem2No) != curvedElements_.end()) {
							throw std::runtime_error("TFSF defined on curved elements is not supported.");
						}
					}
				}
			}
		}
	}
}

bool HesthavenEvolution::isDoFinCurvedElement(const NodeId& d) const
{
	for (int c = 0; c < hestElemCurvedStorage_.size(); c++)
	{
		if (std::find(hestElemCurvedStorage_[c].dofs.begin(), hestElemCurvedStorage_[c].dofs.end(), d) != hestElemCurvedStorage_[c].dofs.end()) {
			return true;
		}
	}
	return false;
}

void HesthavenEvolution::applyBoundaryConditionsToNodes(const BoundaryMaps& bdrMaps,
                                                        const FieldsInputMaps& in,
                                                        HesthavenFields& out) const
{
	auto applyBdrBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < bc.vmapB.size(); m++) {
			for (int d = X; d <= Z; d++) {
				for (int v = 0; v < bc.vmapB[m].size(); v++) {
					if (!isDoFinCurvedElement(bc.vmapB[m][v])) {
						out.e_[d][bc.mapB[m][v]] = e_coeff * in.e_[d][bc.vmapB[m][v]];
						out.h_[d][bc.mapB[m][v]] = h_coeff * in.h_[d][bc.vmapB[m][v]];
					}
				}
			}
		}
	};

	applyBdrBC(bdrMaps.PEC, -2.0, 0.0);
	applyBdrBC(bdrMaps.PMC, 0.0, -2.0);
	if (opts_.alpha != 0.0) {
		applyBdrBC(bdrMaps.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);
	}

	auto applyIntBdrBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < bc.mapBElem1.size(); m++) {
			for (int d = X; d <= Z; d++) {
				for (int v = 0; v < bc.mapBElem1[m].size(); v++) {
					if (!isDoFinCurvedElement(bc.vmapBElem1[m][v]) || !isDoFinCurvedElement(bc.vmapBElem2[m][v])) {
						out.e_[d][bc.mapBElem1[m][v]] = e_coeff * in.e_[d][bc.vmapBElem1[m][v]];
						out.h_[d][bc.mapBElem1[m][v]] = h_coeff * in.h_[d][bc.vmapBElem1[m][v]];
						out.e_[d][bc.mapBElem2[m][v]] = e_coeff * in.e_[d][bc.vmapBElem2[m][v]];
						out.h_[d][bc.mapBElem2[m][v]] = h_coeff * in.h_[d][bc.vmapBElem2[m][v]];
					}
				}
			}
		}
	};

	applyIntBdrBC(bdrMaps.intPEC, -1.0, 0.0);
	applyIntBdrBC(bdrMaps.intPMC, 0.0, -1.0);
	applyIntBdrBC(bdrMaps.intSMA, -0.5, -0.5);
}

void HesthavenEvolution::applyBoundaryConditionsToNodes(const BoundaryMaps& bdrMaps,
                                                        const FieldsInputMaps& in,
                                                        double* e_jumps,
                                                        double* h_jumps,
                                                        int jumps_size) const
{
	auto applyBdrBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < bc.vmapB.size(); m++) {
			for (int d = X; d <= Z; d++) {
				for (int v = 0; v < bc.vmapB[m].size(); v++) {
					if (!isDoFinCurvedElement(bc.vmapB[m][v])) {
						e_jumps[d * jumps_size + bc.mapB[m][v]] = e_coeff * in.e_[d][bc.vmapB[m][v]];
						h_jumps[d * jumps_size + bc.mapB[m][v]] = h_coeff * in.h_[d][bc.vmapB[m][v]];
					}
				}
			}
		}
	};

	applyBdrBC(bdrMaps.PEC, -2.0, 0.0);
	applyBdrBC(bdrMaps.PMC, 0.0, -2.0);
	if (opts_.alpha != 0.0) {
		applyBdrBC(bdrMaps.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);
	}

	auto applyIntBdrBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < bc.mapBElem1.size(); m++) {
			for (int d = X; d <= Z; d++) {
				for (int v = 0; v < bc.mapBElem1[m].size(); v++) {
					if (!isDoFinCurvedElement(bc.vmapBElem1[m][v]) || !isDoFinCurvedElement(bc.vmapBElem2[m][v])) {
						e_jumps[d * jumps_size + bc.mapBElem1[m][v]] = e_coeff * in.e_[d][bc.vmapBElem1[m][v]];
						h_jumps[d * jumps_size + bc.mapBElem1[m][v]] = h_coeff * in.h_[d][bc.vmapBElem1[m][v]];
						e_jumps[d * jumps_size + bc.mapBElem2[m][v]] = e_coeff * in.e_[d][bc.vmapBElem2[m][v]];
						h_jumps[d * jumps_size + bc.mapBElem2[m][v]] = h_coeff * in.h_[d][bc.vmapBElem2[m][v]];
					}
				}
			}
		}
	};

	applyIntBdrBC(bdrMaps.intPEC, -1.0, 0.0);
	applyIntBdrBC(bdrMaps.intPMC, 0.0, -1.0);
	applyIntBdrBC(bdrMaps.intSMA, -0.5, -0.5);

}

HesthavenEvolution::~HesthavenEvolution() = default;

void HesthavenEvolution::SetTime(double t)
{
	TimeDependentOperator::SetTime(t);
	if (sgbc_engine_) {
		sgbc_engine_->SetTime(t);
	}
}

void HesthavenEvolution::applySGBC(const Vector& in, Vector& out, bool inject) const
{
	if (sgbc_engine_ == nullptr || !sgbc_engine_->hasSGBC()) {
		return;
	}
	sgbc_engine_->SetTime(GetTime());
	if (inject) {
		sgbc_engine_->injectSGBCFlux(in, out);
	} else {
		sgbc_engine_->advanceSGBC(in);
	}
}

void HesthavenEvolution::promoteAssembledElements()
{
	auto marked = [](const Array<int>& marker, int attr) {
		return attr > 0 && marker.Size() >= attr && marker[attr - 1] == 1;
	};
	std::unordered_set<int> extra;
	auto add_boundary_element = [&](int boundary) {
		const int face = fes_.GetMesh()->GetBdrElementFaceIndex(boundary);
		const bool interior = fes_.GetMesh()->FaceIsInterior(face);
		const mfem::FaceElementTransformations* ft = interior
			? fes_.GetMesh()->GetInternalBdrFaceTransformations(boundary)
			: fes_.GetMesh()->GetBdrFaceTransformations(boundary);
		if (ft == nullptr) {
			return;
		}
		extra.insert(ft->Elem1No);
		if (interior && ft->Elem2No >= 0 && ft->Elem2No < fes_.GetMesh()->GetNE()) {
			extra.insert(ft->Elem2No);
		}
	};

	const int nbe = fes_.GetMesh()->GetNBE();
	std::unordered_set<int> tfsf_elements;
	for (const auto& entry : model_.getTotalFieldScatteredFieldToMarker()) {
		for (int b = 0; b < nbe; ++b) {
			if (!marked(entry.second, fes_.GetMesh()->GetBdrAttribute(b))) {
				continue;
			}
			const int face = fes_.GetMesh()->GetBdrElementFaceIndex(b);
			const bool interior = fes_.GetMesh()->FaceIsInterior(face);
			const mfem::FaceElementTransformations* ft = interior
				? fes_.GetMesh()->GetInternalBdrFaceTransformations(b)
				: fes_.GetMesh()->GetBdrFaceTransformations(b);
			if (ft == nullptr) {
				continue;
			}
			tfsf_elements.insert(ft->Elem1No);
			if (interior && ft->Elem2No >= 0 && ft->Elem2No < fes_.GetMesh()->GetNE()) {
				tfsf_elements.insert(ft->Elem2No);
			}
		}
	}
	auto blocked_by_tfsf = [&](int boundary) {
		const int face = fes_.GetMesh()->GetBdrElementFaceIndex(boundary);
		const bool interior = fes_.GetMesh()->FaceIsInterior(face);
		const mfem::FaceElementTransformations* ft = interior
			? fes_.GetMesh()->GetInternalBdrFaceTransformations(boundary)
			: fes_.GetMesh()->GetBdrFaceTransformations(boundary);
		if (ft == nullptr) {
			return true;
		}
		if (tfsf_elements.count(ft->Elem1No) != 0) {
			return true;
		}
		return interior && ft->Elem2No >= 0 && tfsf_elements.count(ft->Elem2No) != 0;
	};

	for (int b = 0; b < nbe; ++b) {
		const int attr = fes_.GetMesh()->GetBdrAttribute(b);
		bool ibfi = false;
		bool sgbc = false;
		for (const auto& entry : model_.getInteriorBoundaryToMarker()) {
			if (!marked(entry.second, attr)) {
				continue;
			}
			if (entry.first == BdrCond::PEC || entry.first == BdrCond::PMC || entry.first == BdrCond::SMA) {
				ibfi = true;
			}
			if (entry.first == BdrCond::SGBC) {
				sgbc = true;
			}
		}
		{
			const auto it = model_.getBoundaryToMarker().find(BdrCond::SGBC);
			if (it != model_.getBoundaryToMarker().end() && marked(it->second, attr)) {
				sgbc = true;
			}
		}
		bool centered_sma = false;
		if (opts_.alpha == 0.0) {
			const auto it = model_.getBoundaryToMarker().find(BdrCond::SMA);
			if (it != model_.getBoundaryToMarker().end() && marked(it->second, attr)) {
				centered_sma = true;
			}
		}
		if ((ibfi || centered_sma || sgbc) && !blocked_by_tfsf(b)) {
			add_boundary_element(b);
		}
	}
	if (extra.empty()) {
		return;
	}

	Array<int> kept;
	for (int i = 0; i < linearElements_.Size(); ++i) {
		if (extra.count(linearElements_[i]) == 0) {
			kept.Append(linearElements_[i]);
		}
	}
	linearElements_ = kept;
	for (int elem : extra) {
		if (curvedElements_.find(elem) != curvedElements_.end()) {
			continue;
		}
		Array<int> dofs;
		fes_.GetElementDofs(elem, dofs);
		for (int i = 0; i < dofs.Size(); ++i) {
			if (dofs[i] < 0) {
				dofs[i] = -1 - dofs[i];
			}
		}
		curvedElements_.emplace(elem, dofs);
	}
}

void HesthavenEvolution::ImplicitSolve(const double dt, const Vector& x, Vector& k)
{
	if (sgbc_engine_ != nullptr && sgbc_engine_->hasSGBC()) {
		throw std::runtime_error(
			"SGBC is integrated inside explicit RK4. Use ode_type RK4.");
	}
	const int n = Height();
	const bool on_device = GetMemoryClass() == mfem::MemoryClass::DEVICE;

	if (implicit_rhs_.Size() != n) {
		implicit_rhs_.SetSize(n);
		implicit_rhs_.UseDevice(on_device);
		implicit_tmp_.SetSize(n);
		implicit_tmp_.UseDevice(on_device);
	}
	if (k.Size() != n) {
		k.SetSize(n);
		k.UseDevice(on_device);
	}

	// Inhomogeneous sources: Mult(0). Reused as RHS bias; homogeneous A via Mult(v)-source.
	implicit_tmp_ = 0.0;
	Mult(implicit_tmp_, implicit_rhs_);
	auto apply_a = [&](const Vector& v, Vector& av) {
		Mult(v, av);
		av -= implicit_rhs_;
	};
	apply_a(x, k);
	k += implicit_rhs_;

	class HomogeneousJacobian final : public mfem::Operator {
		std::function<void(const Vector&, Vector&)> apply_;
		const double dt_;
		mfem::Vector& tmp_;
	public:
		HomogeneousJacobian(int size, std::function<void(const Vector&, Vector&)> apply,
		                    double step, mfem::Vector& tmp)
			: mfem::Operator(size), apply_(std::move(apply)), dt_(step), tmp_(tmp) {}
		void Mult(const Vector& v, Vector& jv) const override
		{
			apply_(v, tmp_);
			jv = v;
			jv.Add(-dt_, tmp_);
		}
	} jacobian(n, apply_a, dt, implicit_tmp_);

	mfem::GMRESSolver gmres;
	gmres.SetOperator(jacobian);
	gmres.SetRelTol(1e-8);
	gmres.SetMaxIter(400);
	gmres.SetKDim(60);
	gmres.SetPrintLevel(0);
	implicit_cached_dt_ = dt;
	Vector rhs(k);
	gmres.Mult(rhs, k);
}

HesthavenEvolution::HesthavenEvolution(ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts,
                                       const Probes* probes, double final_time) :
	TimeDependentOperator(hesthavenStateSize(fes, model)),
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts)
{

	const auto* cmesh = &model_.getConstMesh();
	ParMesh& mesh = model_.getMesh();
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };

	mfem::StopWatch swPositions, swClassify, swRefLift, swElemLoop, swDirCrit,
		swGroups, swPromoted, swConnectivity, swTfsfPre, swAux, swSgbc, swGpuInit;
	std::atomic<double> face_info_secs{0.0};

	swPositions.Start();
	positions_ = buildDoFPositions(fes_);
	swPositions.Stop();
	swClassify.Start();
	auto elemOrderList = initCurvedAndLinearElementsLists(fes_, positions_);
	linearElements_ = elemOrderList.first;
	curvedElements_ = elemOrderList.second;
	promoteAssembledElements();
	checkForTFSFInCurvedElements();
	swClassify.Stop();

	hestElemLinearStorage_.resize(linearElements_.Size());

	bool allElementsSameGeomType = true;
	{
		const auto firstElemGeomType = cmesh->GetElementGeometry(0);
		for (auto e = 0; e < cmesh->GetNE(); e++) {
			if (firstElemGeomType != cmesh->GetElementGeometry(e)) {
				allElementsSameGeomType = false;
			}
		}
	}

	swRefLift.Start();
	DynamicMatrix refInvMass;
	if (allElementsSameGeomType) {
		refInvMass = assembleHesthavenRefElemInvMassMatrix(
			cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
		refLIFT_ = refInvMass * assembleHesthavenRefElemEmat(
			cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
	}
	swRefLift.Stop();

	const int n_assemble = linearElements_.Size();
	if (Mpi::WorldRank() == 0) {
		std::cout << "Hesthaven init: assembling " << n_assemble
		          << " linear elements (parent-mesh affine)...\n"
		          << std::flush;
	}

	std::atomic<int> assembled_count{0};
	const int progress_stride = std::max(1, n_assemble / 16);

	auto assembleOneElement = [&](int e) {
		HesthavenElement hestElem;
		hestElem.id = linearElements_[e];
		hestElem.type = cmesh->GetElementType(linearElements_[e]);
		hestElem.vol = mesh.GetElementVolume(linearElements_[e]);

#ifdef SEMBA_DGTD_ENABLE_OPENMP
		#pragma omp critical(hesthaven_matrix_dedup)
#endif
		{
			swDirCrit.Start();
			storeDirectionalMatrices(refInvMass, hestElem);
			swDirCrit.Stop();
		}
		{
			mfem::StopWatch swFace;
			swFace.Start();
			storeFaceInformation(hestElem);
			swFace.Stop();
			face_info_secs.store(face_info_secs.load(std::memory_order_relaxed) + swFace.RealTime(),
			                     std::memory_order_relaxed);
		}
		fes_.GetElementDofs(hestElem.id, hestElem.dofs);
		for (int i = 0; i < hestElem.dofs.Size(); ++i) {
			if (hestElem.dofs[i] < 0) {
				hestElem.dofs[i] = -1 - hestElem.dofs[i];
			}
		}
		hestElemLinearStorage_[e] = std::move(hestElem);

		const int done = ++assembled_count;
		if (Mpi::WorldRank() == 0 && (done % progress_stride == 0 || done == n_assemble)) {
			std::cout << "Hesthaven init: " << done << "/" << n_assemble << " elements\n"
			          << std::flush;
		}
	};

	// The per-element assembly touches MFEM mesh/FES accessors with lazy caches
	// (e.g. face transformations) that are not thread-safe; keep it serial. It
	// measures <1 s for ~59k elements, so parallelizing it is not worth the race.
	swElemLoop.Start();
	for (int e = 0; e < linearElements_.Size(); e++) {
		assembleOneElement(e);
	}
	swElemLoop.Stop();

	swGroups.Start();
	buildLinearDirGroups();
	swGroups.Stop();

	// buildGlobalOperator() starts with ExchangeFaceNbrData — every rank must enter
	// when any rank has promoted/curved elements (SGBC / IBFI / centered SMA), or MPI deadlocks.
	swPromoted.Start();
	{
		const MPI_Comm comm = fes_.GetParMesh()->GetComm();
		int need_assembled = curvedElements_.empty() ? 0 : 1;
		if (Mpi::WorldSize() > 1) {
			MPI_Allreduce(MPI_IN_PLACE, &need_assembled, 1, MPI_INT, MPI_MAX, comm);
		}
		if (need_assembled) {
			Probes probes;
			ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
			DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes_);
			auto global = dgops.buildGlobalOperator();

			const int ndofs = fes_.GetNDofs();
			const int nbrDofs = fes_.num_face_nbr_dofs;
			const int blockSize = ndofs + nbrDofs;
			const int localSize = numberOfFieldComponents * numberOfMaxDimensions * ndofs;
			const int fullSize = numberOfFieldComponents * numberOfMaxDimensions * blockSize;

			for (const auto& [e, dofs] : curvedElements_) {
				HesthavenCurvedElement hestCurElem;
				hestCurElem.id = e;
				hestCurElem.type = cmesh->GetElementType(e);
				hestCurElem.dofs = dofs;
				SparseMatrix spmat(localSize, fullSize);
				Array<int> cols;
				Vector vals;
				for (auto d : dofs) {
					for (auto ft = 0; ft < numberOfFieldComponents * numberOfMaxDimensions; ft++) {
						global->GetRow(d + ft * ndofs, cols, vals);
						spmat.SetRow(d + ft * ndofs, cols, vals);
					}
				}
				spmat.Finalize();
				hestCurElem.matrix = std::move(spmat);
				hestElemCurvedStorage_.push_back(hestCurElem);
			}
		}
	}
	swPromoted.Stop();

	{
		const int P = Mpi::WorldSize();
		const MPI_Comm comm = fes_.GetParMesh()->GetComm();
		const int local_linear = linearElements_.Size();
		const int local_curved = static_cast<int>(curvedElements_.size());
		const int local_unique = static_cast<int>(matrixStorage_.size());

		int global_linear = local_linear;
		int global_curved = local_curved;
		int max_unique = local_unique;
		int min_unique = local_unique;
		if (P > 1) {
			MPI_Allreduce(MPI_IN_PLACE, &global_linear, 1, MPI_INT, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, &global_curved, 1, MPI_INT, MPI_SUM, comm);
			MPI_Allreduce(&local_unique, &max_unique, 1, MPI_INT, MPI_MAX, comm);
			MPI_Allreduce(&local_unique, &min_unique, 1, MPI_INT, MPI_MIN, comm);
		}

		if (Mpi::WorldRank() == 0) {
			const int dir_refs = global_linear * 3;
			const double savings_pct = (dir_refs > 0)
				? 100.0 * (1.0 - static_cast<double>(max_unique) / dir_refs) : 0.0;
			std::cout << "\nHesthaven matrix reuse: " << max_unique
			          << " unique directional matrices for " << global_linear
			          << " linear elements (" << dir_refs << " references, "
			          << std::fixed << std::setprecision(1) << savings_pct << "% saved)";
			if (global_curved > 0) {
				std::cout << ", " << global_curved << " curved (no sharing)";
			}
			if (P > 1 && min_unique != max_unique) {
				std::cout << " [unique/rank: " << min_unique << "-" << max_unique << "]";
			}
			std::cout << "\n" << std::flush;
		}
	}

	fes_.ExchangeFaceNbrData();
	for (int d = X; d <= Z; d++) {
		eOld_[d].SetSpace(&fes_);
		hOld_[d].SetSpace(&fes_);
	}

	swConnectivity.Start();
	connectivity_.emplace(model, fes_);
	swConnectivity.Stop();

	{
		const int n_threads =
#ifdef SEMBA_DGTD_ENABLE_OPENMP
			omp_get_max_threads();
#else
			1;
#endif
		mult_scratch_.resize(static_cast<size_t>(n_threads));
		mult_jumps_ = std::make_unique<HesthavenFields>(connectivity_->global.size());
	}

	swTfsfPre.Start();
	const auto& tfsf_markers = model_.getTotalFieldScatteredFieldToMarker();
	if (!tfsf_markers.empty()) {
		const auto it = tfsf_markers.find(BdrCond::TotalFieldIn);
		if (it != tfsf_markers.end()) {
			srcmngr_.initTFSFPreReqs(model_.getConstMesh(), it->second);
			srcmngr_.initDirectPlanewaveEval();
#ifdef SEMBA_DGTD_ENABLE_CUDA
			has_tfsf_gpu_ = srcmngr_.hasDirectEval();
#endif
		}
	}
	swTfsfPre.Stop();

	swAux.Start();
	auxiliary_.build(fes_, model_, linearElements_);
	swAux.Stop();
	swSgbc.Start();
	if (!model_.getSGBCProperties().empty()) {
		static const Probes empty_probes;
		const Probes& sgbc_probes = probes != nullptr ? *probes : empty_probes;
		sgbc_host_opts_ = std::make_unique<EvolutionOptions>(opts_);
		sgbc_host_opts_->hesthaven_sgbc_host = true;
		sgbc_engine_ = std::make_unique<GlobalEvolution>(
			fes_, model_, srcmngr_, *sgbc_host_opts_, sgbc_probes, final_time);
	}
	swSgbc.Stop();

	buildMergedCurvedMatrix();
#ifdef SEMBA_DGTD_ENABLE_CUDA
	swGpuInit.Start();
	if (linearElements_.Size() > 0 && refLIFT_.size() > 0) {
		initGPUData();
	}
	if (!hestElemCurvedStorage_.empty()) {
		initGPUCurvedData();
	}
	swGpuInit.Stop();
#endif

#ifdef SHOW_TIMER_INFORMATION
	{
		if (Mpi::WorldRank() == 0) {
			std::cout << "------------------------------------------------\n"
			          << "Hesthaven init phases (rank 0, seconds):\n"
			          << "  positions/classify: " << std::fixed << std::setprecision(2)
			          << swPositions.RealTime() + swClassify.RealTime() << "\n"
			          << "  ref inv-mass/lift : " << swRefLift.RealTime() << "\n"
			          << "  element loop      : " << swElemLoop.RealTime()
			          << "  (dir assy+dedup " << swDirCrit.RealTime()
			          << ", face info " << face_info_secs.load() << ")\n"
			          << "  dir groups        : " << swGroups.RealTime() << "\n"
			          << "  promoted/curved   : " << swPromoted.RealTime() << "\n"
			          << "  connectivity      : " << swConnectivity.RealTime() << "\n"
			          << "  tfsf pre-reqs     : " << swTfsfPre.RealTime() << "\n"
			          << "  auxiliary build   : " << swAux.RealTime() << "\n"
			          << "  sgbc host engine  : " << swSgbc.RealTime() << "\n";
#ifdef SEMBA_DGTD_ENABLE_CUDA
			std::cout << "  gpu init          : " << swGpuInit.RealTime() << "\n";
#endif
			std::cout << "------------------------------------------------\n" << std::flush;
		}
	}
	{
		const int P = Mpi::WorldSize();
		const MPI_Comm comm = fes_.GetParMesh()->GetComm();
		const int local_counts[3] = {
			linearElements_.Size(),
			static_cast<int>(curvedElements_.size()),
			cmesh->GetNE()
		};

		std::vector<int> all_counts(3 * P);
		if (P > 1) {
			MPI_Gather(local_counts, 3, MPI_INT, all_counts.data(), 3, MPI_INT, 0, comm);
		} else {
			std::copy(local_counts, local_counts + 3, all_counts.data());
		}

		int global_totals[3] = {local_counts[0], local_counts[1], local_counts[2]};
		if (P > 1) {
			MPI_Allreduce(MPI_IN_PLACE, global_totals, 3, MPI_INT, MPI_SUM, comm);
		}

		if (Mpi::WorldRank() == 0) {
			std::cout << "------------------------------------------------\n\n"
			          << "Hesthaven Evolution Operator initialized.\n\n"
			          << "  Rank  | linear  curved | local total |   load\n"
			          << "  ------+----------------------------+--------\n";
			for (int r = 0; r < P; ++r) {
				const int n_lin = all_counts[r * 3 + 0];
				const int n_cur = all_counts[r * 3 + 1];
				const int n_loc = all_counts[r * 3 + 2];
				const double load_pct = global_totals[2] > 0
					? 100.0 * n_loc / global_totals[2] : 0.0;

				char buf[256];
				std::snprintf(buf, sizeof(buf),
				              "  %4d  | %6d  %6d | %11d | %5.1f%%",
				              r, n_lin, n_cur, n_loc, load_pct);
				std::cout << buf << '\n';
			}
			std::cout << "  ------+----------------------------+--------\n";
			char total_buf[256];
			std::snprintf(total_buf, sizeof(total_buf),
			              "  total | %6d  %6d | %11d |\n\n",
			              global_totals[0], global_totals[1], global_totals[2]);
			std::cout << total_buf
			          << "------------------------------------------------\n" << std::flush;
		}
	}
#endif

}

void loadOutVectors(const Eigen::VectorXd& data, const FiniteElementSpace& fes, const ElementId& e, GridFunction& out)
{
	Array<int> dofs;
	auto el2dofs = fes.GetElementDofs(e, dofs);
	std::unique_ptr<mfem::real_t[]> mfemFieldVars = std::make_unique<mfem::real_t[]>(data.size());
	for (int v = 0; v < data.size(); v++) {
		mfemFieldVars.get()[v] = data.data()[v];
	}
	out.SetSubVector(dofs, mfemFieldVars.get());
}

void HesthavenEvolution::exchangeFieldData(const Vector& in, bool for_gpu_mult) const
{
	const int ndofs = fes_.GetNDofs();
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (for_gpu_mult && gpu_.initialized && hesthavenGpuMultViable()) {
		load_in_to_eh_gpu(in, eOld_, hOld_, ndofs);
	} else
#endif
	{
		const double* in_data = in.HostRead();
		for (int d = X; d <= Z; ++d) {
			double* e_data = eOld_[d].HostWrite();
			double* h_data = hOld_[d].HostWrite();
			std::memcpy(e_data, in_data + d * ndofs, ndofs * sizeof(double));
			std::memcpy(h_data, in_data + (3 + d) * ndofs, ndofs * sizeof(double));
		}
#ifdef SEMBA_DGTD_ENABLE_CUDA
		if (mfem::Device::Allows(mfem::Backend::CUDA)) {
			// HostWrite above; push to device (Write() would not copy and MPI send
			// in ExchangeFaceNbrData reads stale device DOFs).
			for (int d = X; d <= Z; ++d) {
				(void)eOld_[d].Read();
				(void)hOld_[d].Read();
			}
			MFEM_STREAM_SYNC;
		}
#endif
	}
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (for_gpu_mult && gpu_.initialized && hesthavenGpuMultViable()) {
		MFEM_STREAM_SYNC;
	}
#endif
	for (int d = X; d <= Z; ++d) {
		eOld_[d].ExchangeFaceNbrData();
		hOld_[d].ExchangeFaceNbrData();
	}
}

void HesthavenEvolution::computeJumps(const Vector& in, HesthavenFields& jumps) const
{
	const FieldsInputMaps fieldsIn(in, fes_);

	const auto& plus_is_nbr = connectivity_->global_plus_is_face_nbr;
	const bool use_nbr_flags = (plus_is_nbr.size() == connectivity_->global.size());

	for (int v = 0; v < connectivity_->global.size(); v++) {
		const bool from_nbr = use_nbr_flags && plus_is_nbr[v];
		for (int d = X; d <= Z; d++) {
			const int minus = connectivity_->global[v].first;
			const int plus = connectivity_->global[v].second;
			const double e_plus = from_nbr
				? eOld_[d].FaceNbrData().HostRead()[plus]
				: fieldsIn.e_[d][plus];
			const double h_plus = from_nbr
				? hOld_[d].FaceNbrData().HostRead()[plus]
				: fieldsIn.h_[d][plus];
			jumps.e_[d][v] = e_plus - fieldsIn.e_[d][minus];
			jumps.h_[d][v] = h_plus - fieldsIn.h_[d][minus];
		}
	}
}

void HesthavenEvolution::Mult(const Vector& in, Vector& out) const
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (gpu_.initialized && hesthavenGpuMultViable()) {
		MultGPU(in, out);
		return;
	}
#endif
	MultCPU(in, out);
}

mfem::MemoryClass HesthavenEvolution::GetMemoryClass() const
{
#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (gpu_.initialized && hesthavenGpuMultViable()) {
		return mfem::MemoryClass::DEVICE;
	}
#endif
	return mfem::MemoryClass::HOST;
}

void HesthavenEvolution::MultCPU(const Vector& in, Vector& out) const
{
#ifdef SHOW_TIMER_INFORMATION
	mfem::StopWatch timerTotal, timerExchange, timerJumps, timerElements;
	mfem::StopWatch timerPromoted, timerAux, timerSgbc;
	timerTotal.Start();
	double tfsf_abs_sum_e = 0.0;
	double tfsf_abs_sum_h = 0.0;
	double tfsf_max_e = 0.0;
	double tfsf_max_h = 0.0;
	double out_abs_sum = 0.0;
	double out_max = 0.0;
#endif

	double alpha = opts_.alpha;
	const double* in_data = in.HostRead();
	applySGBC(in, out, false);

	out.SetSize(Height());

#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Start();
#endif
	exchangeFieldData(in, false);
#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Stop();
	timerJumps.Start();
#endif

	const FieldsInputMaps fieldsIn(in, fes_);
	FieldsOutputMaps outMaps(out, fes_);

	const int ndofs = fes_.GetNDofs();
	for (int d = X; d <= Z; ++d) {
		outMaps.e_[d].setZero();
		outMaps.h_[d].setZero();
	}
	if (Height() > 6 * ndofs) {
		double* out_tail = out.HostReadWrite();
		std::fill(out_tail + 6 * ndofs, out_tail + static_cast<int>(Height()), 0.0);
	}

	auto scatterElem = [&](Eigen::Ref<const Eigen::VectorXd> data, const Array<int>& dofs,
	                       Eigen::Map<Eigen::VectorXd>& out_map) {
		for (int i = 0; i < dofs.Size(); ++i) {
			out_map[dofs[i]] += data[i];
		}
	};

	auto& jumps = *mult_jumps_;
	computeJumps(in, jumps);

	applyBoundaryConditionsToNodes(connectivity_->boundary, fieldsIn, jumps);
	evaluateTFSF(jumps);
#ifdef SHOW_TIMER_INFORMATION
	for (int d = X; d <= Z; ++d) {
		for (int i = 0; i < jumps.e_[d].size(); ++i) {
			const double ae = std::abs(jumps.e_[d][i]);
			const double ah = std::abs(jumps.h_[d][i]);
			tfsf_abs_sum_e += ae;
			tfsf_abs_sum_h += ah;
			tfsf_max_e = std::max(tfsf_max_e, ae);
			tfsf_max_h = std::max(tfsf_max_h, ah);
		}
	}
#endif
#ifdef SHOW_TIMER_INFORMATION
	timerJumps.Stop();
	timerElements.Start();
#endif

	auto processOneLinear = [&](int e) {
		const auto& he = hestElemLinearStorage_[static_cast<size_t>(e)];
		const int nd = he.dofs.Size();
		const int nf = static_cast<int>(he.fscale.size());
		int tid = 0;
#ifdef SEMBA_DGTD_ENABLE_OPENMP
		tid = omp_get_thread_num();
#endif
		auto& sc = mult_scratch_[static_cast<size_t>(tid)];
		sc.ensure(nd, nf);

		// Gather this element's field values (replaces FieldsElementMaps copies).
		for (int i = 0; i < nd; ++i) {
			const int dof = he.dofs[i];
			for (int d = X; d <= Z; ++d) {
				sc.eIn[d][i] = in_data[d * ndofs + dof];
				sc.hIn[d][i] = in_data[(3 + d) * ndofs + dof];
			}
		}

		// Jump views into the preallocated element-major jump buffer.
		const int jbase = he.id * nf;

		sc.ndotdH.setZero();
		sc.ndotdE.setZero();
		for (int d = X; d <= Z; ++d) {
			const Eigen::Map<const Eigen::VectorXd> jh_d(jumps.h_[d].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> je_d(jumps.e_[d].data() + jbase, nf);
			sc.ndotdH.noalias() += he.normals[d].asDiagonal() * jh_d;
			sc.ndotdE.noalias() += he.normals[d].asDiagonal() * je_d;
		}

		const HesthavenLinearMaterial& mat = auxiliary_.linearMaterial(e);

		for (int x = X; x <= Z; ++x) {
			const int y = (x + 1) % 3;
			const int z = (x + 2) % 3;

			const Eigen::VectorXd& norx = he.normals[x];
			const Eigen::VectorXd& nory = he.normals[y];
			const Eigen::VectorXd& norz = he.normals[z];
			const DynamicMatrix& dir1 = *he.dir[y];
			const DynamicMatrix& dir2 = *he.dir[z];

			const Eigen::Map<const Eigen::VectorXd> je_x(jumps.e_[x].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> jh_x(jumps.h_[x].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> je_y(jumps.e_[y].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> jh_y(jumps.h_[y].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> je_z(jumps.e_[z].data() + jbase, nf);
			const Eigen::Map<const Eigen::VectorXd> jh_z(jumps.h_[z].data() + jbase, nf);

			// Flux h: -nory.*je_z + norz.*je_y + alpha*(jh_x - ndotdH.*norx)
			sc.fluxH.noalias() = nory.asDiagonal() * je_z;
			sc.fluxH *= -1.0;
			sc.tmpF.noalias() = norz.asDiagonal() * je_y;
			sc.fluxH += sc.tmpF;
			sc.tmpF.noalias() = sc.ndotdH.asDiagonal() * norx;
			sc.tmpF2.noalias() = jh_x - sc.tmpF;
			sc.tmpF2 *= alpha;
			sc.fluxH += sc.tmpF2;

			// Flux e: nory.*jh_z - norz.*jh_y + alpha*(je_x - ndotdE.*norx)
			sc.fluxE.noalias() = nory.asDiagonal() * jh_z;
			sc.tmpF.noalias() = norz.asDiagonal() * jh_y;
			sc.tmpF *= -1.0;
			sc.fluxE += sc.tmpF;
			sc.tmpF.noalias() = sc.ndotdE.asDiagonal() * norx;
			sc.tmpF2.noalias() = je_x - sc.tmpF;
			sc.tmpF2 *= alpha;
			sc.fluxE += sc.tmpF2;

			// LIFT (same arithmetic as applyLIFT, in place).
			for (int i = 0; i < nf; ++i) {
				sc.fluxH[i] *= he.fscale[i] / 2.0;
				sc.fluxE[i] *= he.fscale[i] / 2.0;
			}
			sc.liftedH.noalias() = refLIFT_ * sc.fluxH;
			sc.liftedE.noalias() = refLIFT_ * sc.fluxE;

			// Volume h: -dir1*e_z + dir2*e_y + liftedH
			sc.volH.noalias() = dir1 * sc.eIn[z];
			sc.volH *= -1.0;
			sc.tmpV.noalias() = dir2 * sc.eIn[y];
			sc.volH += sc.tmpV;
			sc.volH += sc.liftedH;

			// Volume e: dir1*h_z - dir2*h_y + liftedE
			sc.volE.noalias() = dir1 * sc.hIn[z];
			sc.tmpV.noalias() = dir2 * sc.hIn[y];
			sc.tmpV *= -1.0;
			sc.volE += sc.tmpV;
			sc.volE += sc.liftedE;

			if (mat.active()) {
				if (mat.inv_mu != 1.0) {
					sc.volH *= mat.inv_mu;
				}
				if (mat.inv_eps != 1.0) {
					sc.volE *= mat.inv_eps;
				}
				if (mat.sigma_over_eps != 0.0) {
					sc.volE.noalias() -= mat.sigma_over_eps * sc.eIn[x];
				}
			}
			scatterElem(sc.volH, he.dofs, outMaps.h_[x]);
			scatterElem(sc.volE, he.dofs, outMaps.e_[x]);
		}
	};

	// Batched volume GEMM for groups sharing dir*; flux/LIFT stay per-element.
	constexpr int kBatchMin = 8;
#ifdef SEMBA_DGTD_ENABLE_OPENMP
	#pragma omp parallel for schedule(dynamic) if(linear_dir_groups_.size() > 1)
#endif
	for (int gi = 0; gi < static_cast<int>(linear_dir_groups_.size()); ++gi) {
		const auto& group = linear_dir_groups_[static_cast<size_t>(gi)];
		if (static_cast<int>(group.size()) < kBatchMin) {
			for (int e : group) {
				processOneLinear(e);
			}
			continue;
		}

		const auto& he0 = hestElemLinearStorage_[static_cast<size_t>(group[0])];
		const int n = he0.dofs.Size();
		const int k = static_cast<int>(group.size());
		const int flux_n = static_cast<int>(he0.fscale.size());

		std::vector<FieldsElementMaps> fields;
		std::vector<HesthavenFields> fluxes;
		fields.reserve(static_cast<size_t>(k));
		fluxes.reserve(static_cast<size_t>(k));
		for (int e : group) {
			const auto& he = hestElemLinearStorage_[static_cast<size_t>(e)];
			Eigen::Index flux_size = he.fscale.size();
			const auto jumpsElem{ HesthavenElementJumps(jumps, he.id, flux_size) };
			fields.emplace_back(in, fes_, he.dofs);

			Eigen::VectorXd ndotdH(jumpsElem.h_[X].size()), ndotdE(jumpsElem.e_[X].size());
			ndotdH.setZero(); ndotdE.setZero();
			for (int d = X; d <= Z; d++) {
				ndotdH += he.normals[d].asDiagonal() * jumpsElem.h_[d];
				ndotdE += he.normals[d].asDiagonal() * jumpsElem.e_[d];
			}

			HesthavenFields elemFlux(flux_size);
			for (int x = X; x <= Z; x++) {
				int y = (x + 1) % 3;
				int z = (x + 2) % 3;
				const Eigen::VectorXd& norx = he.normals[x];
				const Eigen::VectorXd& nory = he.normals[y];
				const Eigen::VectorXd& norz = he.normals[z];
				elemFlux.h_[x] = -1.0 * nory.asDiagonal() * jumpsElem.e_[z] +       norz.asDiagonal() * jumpsElem.e_[y] + alpha * (jumpsElem.h_[x] - ndotdH.asDiagonal() * norx);
				elemFlux.e_[x] =        nory.asDiagonal() * jumpsElem.h_[z] - 1.0 * norz.asDiagonal() * jumpsElem.h_[y] + alpha * (jumpsElem.e_[x] - ndotdE.asDiagonal() * norx);
			}
			fluxes.push_back(std::move(elemFlux));
			(void)flux_n;
		}

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;
			const DynamicMatrix& dir1 = *he0.dir[y];
			const DynamicMatrix& dir2 = *he0.dir[z];

			Eigen::MatrixXd Ez(n, k), Ey(n, k), Hz(n, k), Hy(n, k);
			for (int j = 0; j < k; ++j) {
				Ez.col(j) = fields[static_cast<size_t>(j)].e_[z];
				Ey.col(j) = fields[static_cast<size_t>(j)].e_[y];
				Hz.col(j) = fields[static_cast<size_t>(j)].h_[z];
				Hy.col(j) = fields[static_cast<size_t>(j)].h_[y];
			}

			Eigen::MatrixXd volH(n, k), volE(n, k);
			volH.noalias() = -dir1 * Ez + dir2 * Ey;
			volE.noalias() =  dir1 * Hz - dir2 * Hy;

			for (int j = 0; j < k; ++j) {
				const int e = group[static_cast<size_t>(j)];
				const auto& he = hestElemLinearStorage_[static_cast<size_t>(e)];
				const HesthavenLinearMaterial& mat = auxiliary_.linearMaterial(e);
				Eigen::VectorXd fluxH = fluxes[static_cast<size_t>(j)].h_[x];
				Eigen::VectorXd fluxE = fluxes[static_cast<size_t>(j)].e_[x];
				Eigen::VectorXd hResult = volH.col(j) + applyLIFT(he.fscale, fluxH);
				Eigen::VectorXd eResult = volE.col(j) + applyLIFT(he.fscale, fluxE);
				if (mat.active()) {
					if (mat.inv_mu != 1.0) {
						hResult *= mat.inv_mu;
					}
					if (mat.inv_eps != 1.0) {
						eResult *= mat.inv_eps;
					}
					if (mat.sigma_over_eps != 0.0) {
						eResult -= mat.sigma_over_eps * fields[static_cast<size_t>(j)].e_[x];
					}
				}
				scatterElem(hResult, he.dofs, outMaps.h_[x]);
				scatterElem(eResult, he.dofs, outMaps.e_[x]);
			}
		}
	}
#ifdef SHOW_TIMER_INFORMATION
	timerElements.Stop();
	timerPromoted.Start();
#endif

	addCurvedElementContributions(out);
#ifdef SHOW_TIMER_INFORMATION
	timerPromoted.Stop();
	timerAux.Start();
#endif
	auxiliary_.apply(in, out);
#ifdef SHOW_TIMER_INFORMATION
	timerAux.Stop();
	timerSgbc.Start();
#endif
	applySGBC(in, out, true);
#ifdef SHOW_TIMER_INFORMATION
	timerSgbc.Stop();
#endif

#ifdef SHOW_TIMER_INFORMATION
	{
		double elem_out_abs_sum = 0.0;
		double elem_out_max = 0.0;
#ifdef SEMBA_DGTD_ENABLE_CUDA
		if (gpu_.d_elem_out_e.Size() > 0) {
			const double* elem_e = gpu_.d_elem_out_e.HostRead();
			const double* elem_h = gpu_.d_elem_out_h.HostRead();
			for (int i = 0; i < gpu_.d_elem_out_e.Size(); ++i) {
				const double ae = std::abs(elem_e[i]);
				const double ah = std::abs(elem_h[i]);
				elem_out_abs_sum += ae + ah;
				elem_out_max = std::max(elem_out_max, std::max(ae, ah));
			}
		}
#endif
		const double* out_data = out.HostRead();
		for (int i = 0; i < out.Size(); ++i) {
			const double a = std::abs(out_data[i]);
			out_abs_sum += a;
			out_max = std::max(out_max, a);
		}
#ifdef SEMBA_DGTD_ENABLE_CUDA
		if (gpu_.has_tfsf && has_tfsf_gpu_) {
			std::cout << "  elem debug | abs(elem_out)=" << std::scientific << std::setprecision(3)
			          << elem_out_abs_sum
			          << " max(elem_out)=" << elem_out_max
			          << '\n' << std::defaultfloat;
		}
#endif
	}
#endif

#ifdef SHOW_TIMER_INFORMATION
	timerTotal.Stop();

	static double acc_total = 0, acc_exchange = 0, acc_jumps = 0, acc_elements = 0;
	static double acc_promoted = 0, acc_aux = 0, acc_sgbc = 0;
	static int acc_count = 0;
	static double lastPrintTime = -1.0;
	constexpr double printInterval = 0.05;
	constexpr int printCadence = 50;

	acc_total    += timerTotal.RealTime()    * 1000.0;
	acc_exchange += timerExchange.RealTime() * 1000.0;
	acc_jumps    += timerJumps.RealTime()    * 1000.0;
	acc_elements += timerElements.RealTime() * 1000.0;
	acc_promoted += timerPromoted.RealTime() * 1000.0;
	acc_aux      += timerAux.RealTime()      * 1000.0;
	acc_sgbc     += timerSgbc.RealTime()     * 1000.0;
	acc_count++;

	const double currentTime = GetTime();
	bool want_print = (acc_count >= printCadence &&
		(lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval));

	const int P = Mpi::WorldSize();
	if (P > 1) {
		int local_want = want_print ? 1 : 0;
		int global_want = 0;
		MPI_Allreduce(&local_want, &global_want, 1, MPI_INT, MPI_MAX,
		              fes_.GetParMesh()->GetComm());
		want_print = (global_want != 0);
	}

	if (want_print)
	{
		double local_avg[7] = {
			acc_total    / acc_count,
			acc_exchange / acc_count,
			acc_jumps    / acc_count,
			acc_elements / acc_count,
			acc_promoted / acc_count,
			acc_aux      / acc_count,
			acc_sgbc     / acc_count
		};

		std::vector<double> all_avg(7 * P);
		if (P > 1) {
			MPI_Gather(local_avg, 7, MPI_DOUBLE,
			           all_avg.data(), 7, MPI_DOUBLE,
			           0, fes_.GetParMesh()->GetComm());
		} else {
			std::copy(local_avg, local_avg + 7, all_avg.data());
		}

		if (Mpi::WorldRank() == 0) {
			std::cout << "\n[Hesthaven Mult timing] t=" << std::fixed << std::setprecision(4)
			          << currentTime << "  (avg of " << acc_count << " calls, ms/call)\n"
			          << std::defaultfloat;
			std::cout << "  Rank  | total  exchange  jumps  elements  promoted    aux   sgbc  | bottleneck\n"
			          << "  ------+-----------------------------------------------------------+-----------\n";
			for (int r = 0; r < P; ++r) {
				const double t_total    = all_avg[r * 7 + 0];
				const double t_exchange = all_avg[r * 7 + 1];
				const double t_jumps    = all_avg[r * 7 + 2];
				const double t_elements = all_avg[r * 7 + 3];
				const double t_promoted = all_avg[r * 7 + 4];
				const double t_aux      = all_avg[r * 7 + 5];
				const double t_sgbc     = all_avg[r * 7 + 6];
				const char* bottleneck = "elements";
				if (t_exchange > 0.7 * t_total) {
					bottleneck = "MPI wait";
				} else if (t_jumps > 0.5 * t_total) {
					bottleneck = "jumps/BC";
				} else if (t_aux > 0.3 * t_total) {
					bottleneck = "aux/PML";
				} else if (t_promoted > 0.3 * t_total) {
					bottleneck = "promoted";
				} else if (t_sgbc > 0.3 * t_total) {
					bottleneck = "sgbc";
				}

				char buf[320];
				std::snprintf(buf, sizeof(buf),
				              "  %4d  | %6.2f  %8.2f  %6.2f  %8.2f  %8.2f  %6.2f  %5.2f  | %s",
				              r, t_total, t_exchange, t_jumps, t_elements, t_promoted, t_aux, t_sgbc,
				              bottleneck);
				std::cout << buf << '\n';
			}
			std::cout << "  TFSF debug | abs(E)=" << std::scientific << std::setprecision(3)
			          << tfsf_abs_sum_e
			          << " abs(H)=" << tfsf_abs_sum_h
			          << " max(E)=" << tfsf_max_e
			          << " max(H)=" << tfsf_max_h
			          << '\n';
			std::cout << "  out  debug | abs(out)=" << out_abs_sum
			          << " max(out)=" << out_max
			          << '\n' << std::defaultfloat;
			std::cout << std::flush;
		}

		acc_total = acc_exchange = acc_jumps = acc_elements = 0;
		acc_promoted = acc_aux = acc_sgbc = 0;
		acc_count = 0;
		lastPrintTime = currentTime;
	}
#endif
}

#ifdef SEMBA_DGTD_ENABLE_CUDA
void HesthavenEvolution::initGPUData()
{
	if (!mfem::Device::Allows(mfem::Backend::CUDA)) {
		return;
	}
	if (Mpi::WorldSize() > 1 && !mfem::Device::GetGPUAwareMPI() && Mpi::WorldRank() == 0) {
		std::cerr << "Hesthaven CUDA: " << Mpi::WorldSize()
		          << " MPI ranks on a non-GPU-aware MPI build — host-staged halos, "
		          << "duplicated GPU memory per rank.\n";
	}

	const int ndofs = fes_.GetNDofs();
	const int n_elem = linearElements_.Size();
	if (n_elem == 0 || refLIFT_.size() == 0) {
		return;
	}

	gpu_.ndofs = ndofs;
	gpu_.n_elements = n_elem;
	gpu_.jumps_size = connectivity_->global.size();
	gpu_.lift_rows = static_cast<int>(refLIFT_.rows());
	gpu_.lift_cols = static_cast<int>(refLIFT_.cols());

	Array<int> sample_dofs;
	fes_.GetElementDofs(hestElemLinearStorage_[0].id, sample_dofs);
	gpu_.ndof_el = sample_dofs.Size();

	gpu_.max_flux_size = 0;
	for (int e = 0; e < n_elem; ++e) {
		gpu_.max_flux_size = std::max(gpu_.max_flux_size,
		                              static_cast<int>(hestElemLinearStorage_[e].fscale.size()));
	}

	std::vector<const DynamicMatrix*> unique_matrices;
	unique_matrices.reserve(matrixStorage_.size());
	for (auto it = matrixStorage_.begin(); it != matrixStorage_.end(); ++it) {
		unique_matrices.push_back(&(*it));
	}
	auto matrix_index = [&](const DynamicMatrix* mat_ptr) -> int {
		const auto it = matrixStorage_.find(*mat_ptr);
		if (it == matrixStorage_.end()) {
			throw std::runtime_error("Hesthaven GPU init: directional matrix not in storage.");
		}
		return static_cast<int>(std::distance(matrixStorage_.begin(), it));
	};
	gpu_.n_matrices = static_cast<int>(unique_matrices.size());

	size_t total_matrix_entries = 0;
	gpu_.d_matrix_offsets.SetSize(gpu_.n_matrices);
	gpu_.d_matrix_rows.SetSize(gpu_.n_matrices);
	gpu_.d_matrix_cols.SetSize(gpu_.n_matrices);
	gpu_.d_matrix_offsets.GetMemory().UseDevice(true);
	gpu_.d_matrix_rows.GetMemory().UseDevice(true);
	gpu_.d_matrix_cols.GetMemory().UseDevice(true);
	int* matrix_offsets_host = gpu_.d_matrix_offsets.HostWrite();
	int* matrix_rows_host = gpu_.d_matrix_rows.HostWrite();
	int* matrix_cols_host = gpu_.d_matrix_cols.HostWrite();
	for (int m = 0; m < gpu_.n_matrices; ++m) {
		matrix_offsets_host[m] = static_cast<int>(total_matrix_entries);
		matrix_rows_host[m] = static_cast<int>(unique_matrices[m]->rows());
		matrix_cols_host[m] = static_cast<int>(unique_matrices[m]->cols());
		total_matrix_entries += static_cast<size_t>(unique_matrices[m]->rows())
		                      * static_cast<size_t>(unique_matrices[m]->cols());
	}

	gpu_.d_matrices.SetSize(static_cast<int>(total_matrix_entries));
	gpu_.d_matrices.UseDevice(true);
	{
		double* dst = gpu_.d_matrices.HostWrite();
		for (int m = 0; m < gpu_.n_matrices; ++m) {
			const auto& mat = *unique_matrices[m];
			std::memcpy(dst + gpu_.d_matrix_offsets[m], mat.data(),
			            static_cast<size_t>(mat.size()) * sizeof(double));
		}
	}

	gpu_.d_ref_lift.SetSize(gpu_.lift_rows * gpu_.lift_cols);
	gpu_.d_ref_lift.UseDevice(true);
	std::memcpy(gpu_.d_ref_lift.HostWrite(), refLIFT_.data(),
	            static_cast<size_t>(refLIFT_.size()) * sizeof(double));

	gpu_.d_elem_ids.SetSize(n_elem);
	gpu_.d_elem_dofs.SetSize(n_elem * gpu_.ndof_el);
	gpu_.d_dir_matrix_id.SetSize(3 * n_elem);
	gpu_.d_flux_size.SetSize(n_elem);
	gpu_.d_elem_ids.GetMemory().UseDevice(true);
	gpu_.d_elem_dofs.GetMemory().UseDevice(true);
	gpu_.d_dir_matrix_id.GetMemory().UseDevice(true);
	gpu_.d_flux_size.GetMemory().UseDevice(true);
	gpu_.d_normals.SetSize(3 * gpu_.max_flux_size * n_elem);
	gpu_.d_fscale.SetSize(gpu_.max_flux_size * n_elem);
	gpu_.d_normals.UseDevice(true);
	gpu_.d_fscale.UseDevice(true);

	int* elem_ids_host = gpu_.d_elem_ids.HostWrite();
	int* elem_dofs_host = gpu_.d_elem_dofs.HostWrite();
	int* dir_matrix_id_host = gpu_.d_dir_matrix_id.HostWrite();
	int* flux_size_host = gpu_.d_flux_size.HostWrite();
	double* nor_host = gpu_.d_normals.HostWrite();
	double* fsc_host = gpu_.d_fscale.HostWrite();

	for (int e = 0; e < n_elem; ++e) {
		const auto& he = hestElemLinearStorage_[e];
		elem_ids_host[e] = he.id;
		flux_size_host[e] = static_cast<int>(he.fscale.size());
		for (int i = 0; i < he.dofs.Size(); ++i) {
			elem_dofs_host[e * gpu_.ndof_el + i] = he.dofs[i];
		}
		for (int d = 0; d < 3; ++d) {
			dir_matrix_id_host[e * 3 + d] = matrix_index(he.dir[d]);
		}
		for (int d = 0; d < 3; ++d) {
			for (int i = 0; i < he.fscale.size(); ++i) {
				nor_host[e * 3 * gpu_.max_flux_size + d * gpu_.max_flux_size + i] = he.normals[d][i];
			}
		}
		for (int i = 0; i < he.fscale.size(); ++i) {
			fsc_host[e * gpu_.max_flux_size + i] = he.fscale[i];
		}
	}

	gpu_.d_jump_minus.SetSize(gpu_.jumps_size);
	gpu_.d_jump_plus.SetSize(gpu_.jumps_size);
	gpu_.d_jump_plus_is_nbr.SetSize(gpu_.jumps_size);
	gpu_.d_jump_minus.GetMemory().UseDevice(true);
	gpu_.d_jump_plus.GetMemory().UseDevice(true);
	gpu_.d_jump_plus_is_nbr.GetMemory().UseDevice(true);
	int* jmp_minus = gpu_.d_jump_minus.HostWrite();
	int* jmp_plus = gpu_.d_jump_plus.HostWrite();
	uint8_t* jmp_nbr = gpu_.d_jump_plus_is_nbr.HostWrite();
	for (int v = 0; v < gpu_.jumps_size; ++v) {
		jmp_minus[v] = connectivity_->global[v].first;
		jmp_plus[v] = connectivity_->global[v].second;
		jmp_nbr[v] = (v < static_cast<int>(connectivity_->global_plus_is_face_nbr.size()))
			? connectivity_->global_plus_is_face_nbr[v] : 0;
	}

	gpu_.d_jumps_e.SetSize(3 * gpu_.jumps_size);
	gpu_.d_jumps_h.SetSize(3 * gpu_.jumps_size);
	gpu_.d_jumps_e.UseDevice(true);
	gpu_.d_jumps_h.UseDevice(true);
	gpu_.d_elem_out_e.SetSize(3 * n_elem * gpu_.ndof_el);
	gpu_.d_elem_out_h.SetSize(3 * n_elem * gpu_.ndof_el);
	gpu_.d_elem_out_e.UseDevice(true);
	gpu_.d_elem_out_h.UseDevice(true);

	initGPUBoundaryData();

	if (gpu_.ndof_el > 64 || gpu_.max_flux_size > 64 || gpu_.lift_rows > 64) {
		throw std::runtime_error(
			"Hesthaven GPU kernel supports ndof_el, max_flux_size, and lift_rows <= 64.");
	}
	gpu_.workspace_stride = 15 * gpu_.max_flux_size + 9 * gpu_.ndof_el;
	gpu_.d_workspace.SetSize(gpu_.workspace_stride * n_elem);
	gpu_.d_workspace.UseDevice(true);
	gpu_.team_size = std::max({gpu_.ndof_el, gpu_.max_flux_size, gpu_.lift_rows});

	packHesthavenAuxGpu(auxiliary_, gpu_);
	hesthaven_sync_gpu_static_data(gpu_);

	gpu_.initialized = true;
}

void HesthavenEvolution::initGPUCurvedData()
{
	if (!mfem::Device::Allows(mfem::Backend::CUDA)) {
		return;
	}
	if (!curved_merged_matrix_ || curved_merged_matrix_->Height() == 0) {
		return;
	}

	const int n_rows = curved_merged_matrix_->Height();
	const int fullSize = curved_merged_matrix_->Width();

	gpu_.d_curved_row_to_out.SetSize(n_rows);
	gpu_.d_curved_row_to_out.GetMemory().UseDevice(true);
	int* row_to_out = gpu_.d_curved_row_to_out.HostWrite();
	for (int row = 0; row < n_rows; ++row) {
		row_to_out[row] = curved_merged_row_to_out_[row];
	}

	gpu_.n_curved_rows = n_rows;
	gpu_.curved_full_size = fullSize;
	gpu_.d_ext_in.SetSize(fullSize);
	gpu_.d_ext_in.UseDevice(true);
	gpu_.d_curved_y.SetSize(n_rows);
	gpu_.d_curved_y.UseDevice(true);
	(void)gpu_.d_curved_row_to_out.Read();

	gpu_.curved_initialized = true;
}

void HesthavenEvolution::initGPUBoundaryData()
{
	const BoundaryMaps& bdr = connectivity_->boundary;
	const int jumps_size = gpu_.jumps_size;

	std::vector<int> true_jump_out, true_vmap_in, true_comp;
	std::vector<double> true_e_coeff, true_h_coeff;

	auto appendTrueBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < static_cast<int>(bc.vmapB.size()); ++m) {
			for (int d = X; d <= Z; ++d) {
				for (int v = 0; v < static_cast<int>(bc.vmapB[m].size()); ++v) {
					if (!isDoFinCurvedElement(bc.vmapB[m][v])) {
						true_jump_out.push_back(bc.mapB[m][v]);
						true_vmap_in.push_back(bc.vmapB[m][v]);
						true_comp.push_back(d);
						true_e_coeff.push_back(e_coeff);
						true_h_coeff.push_back(h_coeff);
					}
				}
			}
		}
	};

	appendTrueBC(bdr.PEC, -2.0, 0.0);
	appendTrueBC(bdr.PMC, 0.0, -2.0);
	if (opts_.alpha != 0.0) {
		appendTrueBC(bdr.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);
	}

	std::vector<int> int_jump_out1, int_jump_out2, int_vmap_in1, int_vmap_in2, int_comp;
	std::vector<double> int_e_coeff, int_h_coeff;

	auto appendIntBC = [&](const auto& bc, double e_coeff, double h_coeff) {
		for (int m = 0; m < static_cast<int>(bc.mapBElem1.size()); ++m) {
			for (int d = X; d <= Z; ++d) {
				for (int v = 0; v < static_cast<int>(bc.mapBElem1[m].size()); ++v) {
					if (!isDoFinCurvedElement(bc.vmapBElem1[m][v])
					    || !isDoFinCurvedElement(bc.vmapBElem2[m][v])) {
						int_jump_out1.push_back(bc.mapBElem1[m][v]);
						int_jump_out2.push_back(bc.mapBElem2[m][v]);
						int_vmap_in1.push_back(bc.vmapBElem1[m][v]);
						int_vmap_in2.push_back(bc.vmapBElem2[m][v]);
						int_comp.push_back(d);
						int_e_coeff.push_back(e_coeff);
						int_h_coeff.push_back(h_coeff);
					}
				}
			}
		}
	};

	appendIntBC(bdr.intPEC, -1.0, 0.0);
	appendIntBC(bdr.intPMC, 0.0, -1.0);
	appendIntBC(bdr.intSMA, -0.5, -0.5);

	gpu_.n_bc_true = static_cast<int>(true_jump_out.size());
	gpu_.d_bc_true_jump_out.SetSize(gpu_.n_bc_true);
	gpu_.d_bc_true_vmap_in.SetSize(gpu_.n_bc_true);
	gpu_.d_bc_true_comp.SetSize(gpu_.n_bc_true);
	gpu_.d_bc_true_e_coeff.SetSize(gpu_.n_bc_true);
	gpu_.d_bc_true_h_coeff.SetSize(gpu_.n_bc_true);
	gpu_.d_bc_true_jump_out.GetMemory().UseDevice(true);
	gpu_.d_bc_true_vmap_in.GetMemory().UseDevice(true);
	gpu_.d_bc_true_comp.GetMemory().UseDevice(true);
	gpu_.d_bc_true_e_coeff.UseDevice(true);
	gpu_.d_bc_true_h_coeff.UseDevice(true);
	int* bc_true_jump_out = gpu_.d_bc_true_jump_out.HostWrite();
	int* bc_true_vmap_in = gpu_.d_bc_true_vmap_in.HostWrite();
	int* bc_true_comp = gpu_.d_bc_true_comp.HostWrite();
	double* bc_true_e_coeff = gpu_.d_bc_true_e_coeff.HostWrite();
	double* bc_true_h_coeff = gpu_.d_bc_true_h_coeff.HostWrite();
	for (int i = 0; i < gpu_.n_bc_true; ++i) {
		bc_true_jump_out[i] = true_jump_out[i];
		bc_true_vmap_in[i] = true_vmap_in[i];
		bc_true_comp[i] = true_comp[i];
		bc_true_e_coeff[i] = true_e_coeff[i];
		bc_true_h_coeff[i] = true_h_coeff[i];
	}

	gpu_.n_bc_int = static_cast<int>(int_jump_out1.size());
	gpu_.d_bc_int_jump_out1.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_jump_out2.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_vmap_in1.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_vmap_in2.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_comp.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_e_coeff.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_h_coeff.SetSize(gpu_.n_bc_int);
	gpu_.d_bc_int_jump_out1.GetMemory().UseDevice(true);
	gpu_.d_bc_int_jump_out2.GetMemory().UseDevice(true);
	gpu_.d_bc_int_vmap_in1.GetMemory().UseDevice(true);
	gpu_.d_bc_int_vmap_in2.GetMemory().UseDevice(true);
	gpu_.d_bc_int_comp.GetMemory().UseDevice(true);
	gpu_.d_bc_int_e_coeff.UseDevice(true);
	gpu_.d_bc_int_h_coeff.UseDevice(true);
	int* bc_int_jump_out1 = gpu_.d_bc_int_jump_out1.HostWrite();
	int* bc_int_jump_out2 = gpu_.d_bc_int_jump_out2.HostWrite();
	int* bc_int_vmap_in1 = gpu_.d_bc_int_vmap_in1.HostWrite();
	int* bc_int_vmap_in2 = gpu_.d_bc_int_vmap_in2.HostWrite();
	int* bc_int_comp = gpu_.d_bc_int_comp.HostWrite();
	double* bc_int_e_coeff = gpu_.d_bc_int_e_coeff.HostWrite();
	double* bc_int_h_coeff = gpu_.d_bc_int_h_coeff.HostWrite();
	for (int i = 0; i < gpu_.n_bc_int; ++i) {
		bc_int_jump_out1[i] = int_jump_out1[i];
		bc_int_jump_out2[i] = int_jump_out2[i];
		bc_int_vmap_in1[i] = int_vmap_in1[i];
		bc_int_vmap_in2[i] = int_vmap_in2[i];
		bc_int_comp[i] = int_comp[i];
		bc_int_e_coeff[i] = int_e_coeff[i];
		bc_int_h_coeff[i] = int_h_coeff[i];
	}

	std::vector<int> tfsf_jump_sf, tfsf_jump_tf, tfsf_src_dof;
	if (has_tfsf_gpu_) {
		const auto& mapBSF = bdr.TFSF.mapBSF;
		const auto& mapBTF = bdr.TFSF.mapBTF;
		const auto& vmapBSF = bdr.TFSF.vmapBSF;
		for (int m = 0; m < static_cast<int>(mapBSF.size()); ++m) {
			for (int v = 0; v < static_cast<int>(mapBSF[m].size()); ++v) {
				const int pos_idx = vmapBSF[m][v];
				const int src_dof = srcmngr_.findTFSFDofAtPosition(positions_[pos_idx]);
				if (src_dof < 0) {
					throw std::runtime_error(
						"Hesthaven GPU TFSF: could not map jump node to TFSF submesh DOF.");
				}
				tfsf_jump_sf.push_back(mapBSF[m][v]);
				tfsf_jump_tf.push_back(mapBTF[m][v]);
				tfsf_src_dof.push_back(src_dof);
			}
		}
	}

	gpu_.has_tfsf = !tfsf_jump_sf.empty();
	gpu_.n_tfsf = static_cast<int>(tfsf_jump_sf.size());
	gpu_.d_tfsf_jump_sf.SetSize(gpu_.n_tfsf);
	gpu_.d_tfsf_jump_tf.SetSize(gpu_.n_tfsf);
	gpu_.d_tfsf_src_dof.SetSize(gpu_.n_tfsf);
	gpu_.d_tfsf_jump_sf.GetMemory().UseDevice(true);
	gpu_.d_tfsf_jump_tf.GetMemory().UseDevice(true);
	gpu_.d_tfsf_src_dof.GetMemory().UseDevice(true);
	int* tfsf_jump_sf_host = gpu_.d_tfsf_jump_sf.HostWrite();
	int* tfsf_jump_tf_host = gpu_.d_tfsf_jump_tf.HostWrite();
	int* tfsf_src_dof_host = gpu_.d_tfsf_src_dof.HostWrite();
	for (int i = 0; i < gpu_.n_tfsf; ++i) {
		tfsf_jump_sf_host[i] = tfsf_jump_sf[i];
		tfsf_jump_tf_host[i] = tfsf_jump_tf[i];
		tfsf_src_dof_host[i] = tfsf_src_dof[i];
	}

	(void)jumps_size;
}

void HesthavenEvolution::MultGPU(const Vector& in, Vector& out) const
{
#ifdef SHOW_TIMER_INFORMATION
	mfem::StopWatch timerTotal, timerExchange, timerJumps, timerElements;
	mfem::StopWatch timerPromoted, timerAux, timerSgbc;
	timerTotal.Start();
	double tfsf_abs_sum_e = 0.0;
	double tfsf_abs_sum_h = 0.0;
	double tfsf_max_e = 0.0;
	double tfsf_max_h = 0.0;
	double out_abs_sum = 0.0;
	double out_max = 0.0;
	double out_abs_sum_stale = 0.0;
	double out_max_stale = 0.0;
	double elem_out_abs_sum = 0.0;
	double elem_out_max = 0.0;
#endif

	const int ndofs = fes_.GetNDofs();
	applySGBC(in, out, false);

	out.SetSize(Height());
	out.UseDevice(true);
	out = 0.0;

#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Start();
#endif
	exchangeFieldData(in, true);
	if (fes_.num_face_nbr_dofs > 0) {
		sync_cuda_face_nbr_halos(eOld_, hOld_);
	}
#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Stop();
	timerJumps.Start();
#endif

	const_cast<mfem::Vector&>(in).Read();
	hesthaven_compute_jumps_gpu(gpu_, in, eOld_, hOld_, ndofs);
	hesthaven_apply_bc_gpu(gpu_, in, ndofs);
	if (gpu_.has_tfsf && has_tfsf_gpu_) {
		srcmngr_.evalTimeVarFieldDirect(GetTime());
		hesthaven_scatter_tfsf_to_jumps_gpu(
			gpu_, srcmngr_.getCachedTFSFFields(), gpu_.jumps_size);
	} else if (!connectivity_->boundary.TFSF.mapBSF.empty()) {
		// Fallback: host TFSF when direct GPU eval is unavailable.
		(void)in.HostRead();
		auto jumps{ HesthavenFields(gpu_.jumps_size) };
		const double* je = gpu_.d_jumps_e.HostRead();
		const double* jh = gpu_.d_jumps_h.HostRead();
		const int js = gpu_.jumps_size;
		for (int d = X; d <= Z; ++d) {
			std::memcpy(jumps.e_[d].data(), je + d * js, static_cast<size_t>(js) * sizeof(double));
			std::memcpy(jumps.h_[d].data(), jh + d * js, static_cast<size_t>(js) * sizeof(double));
		}
		evaluateTFSF(jumps);
		double* jumps_e_host = gpu_.d_jumps_e.HostWrite();
		double* jumps_h_host = gpu_.d_jumps_h.HostWrite();
		for (int d = X; d <= Z; ++d) {
			std::memcpy(jumps_e_host + d * js, jumps.e_[d].data(),
			            static_cast<size_t>(js) * sizeof(double));
			std::memcpy(jumps_h_host + d * js, jumps.h_[d].data(),
			            static_cast<size_t>(js) * sizeof(double));
		}
		(void)gpu_.d_jumps_e.Read();
		(void)gpu_.d_jumps_h.Read();
	}

#ifdef SHOW_TIMER_INFORMATION
	timerJumps.Stop();
	timerElements.Start();
#endif

	if (linearElements_.Size() > 0) {
		hesthaven_mult_gpu(gpu_, opts_.alpha, in, eOld_, hOld_, out, ndofs);
	}

#ifdef SHOW_TIMER_INFORMATION
	timerElements.Stop();
	timerPromoted.Start();
#endif

	addCurvedElementContributions(out);
#ifdef SHOW_TIMER_INFORMATION
	timerPromoted.Stop();
	timerAux.Start();
#endif

	if (auxiliary_.hasAux()) {
		hesthaven_apply_aux_gpu(gpu_, in, out, ndofs);
	}
#ifdef SHOW_TIMER_INFORMATION
	timerAux.Stop();
	timerSgbc.Start();
#endif
	if (sgbc_engine_ != nullptr && sgbc_engine_->hasSGBC()) {
		(void)out.HostReadWrite();
		applySGBC(in, out, true);
		(void)out.Read();
	}
#ifdef SHOW_TIMER_INFORMATION
	timerSgbc.Stop();
#endif

#ifdef SHOW_TIMER_INFORMATION
	{
		const double* out_data = out.HostRead();
		for (int i = 0; i < out.Size(); ++i) {
			const double a = std::abs(out_data[i]);
			out_abs_sum += a;
			out_max = std::max(out_max, a);
		}
	}
#endif

#ifdef SHOW_TIMER_INFORMATION
	timerTotal.Stop();

	static double acc_total = 0, acc_exchange = 0, acc_jumps = 0, acc_elements = 0;
	static double acc_promoted = 0, acc_aux = 0, acc_sgbc = 0;
	static int acc_count = 0;
	static double lastPrintTime = -1.0;
	constexpr double printInterval = 0.05;
	constexpr int printCadence = 50;

	acc_total    += timerTotal.RealTime()    * 1000.0;
	acc_exchange += timerExchange.RealTime() * 1000.0;
	acc_jumps    += timerJumps.RealTime()    * 1000.0;
	acc_elements += timerElements.RealTime() * 1000.0;
	acc_promoted += timerPromoted.RealTime() * 1000.0;
	acc_aux      += timerAux.RealTime()      * 1000.0;
	acc_sgbc     += timerSgbc.RealTime()     * 1000.0;
	acc_count++;

	const double currentTime = GetTime();
	bool want_print = (acc_count >= printCadence &&
		(lastPrintTime < 0.0 || currentTime >= lastPrintTime + printInterval));

	const int P = Mpi::WorldSize();
	if (P > 1) {
		int local_want = want_print ? 1 : 0;
		int global_want = 0;
		MPI_Allreduce(&local_want, &global_want, 1, MPI_INT, MPI_MAX,
		              fes_.GetParMesh()->GetComm());
		want_print = (global_want != 0);
	}

	if (want_print)
	{
		double local_avg[7] = {
			acc_total    / acc_count,
			acc_exchange / acc_count,
			acc_jumps    / acc_count,
			acc_elements / acc_count,
			acc_promoted / acc_count,
			acc_aux      / acc_count,
			acc_sgbc     / acc_count
		};

		std::vector<double> all_avg(7 * P);
		if (P > 1) {
			MPI_Gather(local_avg, 7, MPI_DOUBLE,
			           all_avg.data(), 7, MPI_DOUBLE,
			           0, fes_.GetParMesh()->GetComm());
		} else {
			std::copy(local_avg, local_avg + 7, all_avg.data());
		}

		if (Mpi::WorldRank() == 0) {
			std::cout << "\n[Hesthaven Mult timing] t=" << std::fixed << std::setprecision(4)
			          << currentTime << "  (avg of " << acc_count << " calls, ms/call, GPU)\n"
			          << std::defaultfloat;
			std::cout << "  Rank  | total  exchange  jumps  elements  promoted    aux   sgbc  | bottleneck\n"
			          << "  ------+-----------------------------------------------------------+-----------\n";
			for (int r = 0; r < P; ++r) {
				const double t_total    = all_avg[r * 7 + 0];
				const double t_exchange = all_avg[r * 7 + 1];
				const double t_jumps    = all_avg[r * 7 + 2];
				const double t_elements = all_avg[r * 7 + 3];
				const double t_promoted = all_avg[r * 7 + 4];
				const double t_aux      = all_avg[r * 7 + 5];
				const double t_sgbc     = all_avg[r * 7 + 6];
				const char* bottleneck = "elements";
				if (t_exchange > 0.7 * t_total) {
					bottleneck = "MPI wait";
				} else if (t_jumps > 0.5 * t_total) {
					bottleneck = "jumps/BC";
				} else if (t_aux > 0.3 * t_total) {
					bottleneck = "aux/PML";
				} else if (t_promoted > 0.3 * t_total) {
					bottleneck = "promoted";
				} else if (t_sgbc > 0.3 * t_total) {
					bottleneck = "sgbc";
				}

				char buf[320];
				std::snprintf(buf, sizeof(buf),
				              "  %4d  | %6.2f  %8.2f  %6.2f  %8.2f  %8.2f  %6.2f  %5.2f  | %s",
				              r, t_total, t_exchange, t_jumps, t_elements, t_promoted, t_aux, t_sgbc,
				              bottleneck);
				std::cout << buf << '\n';
			}
			if (gpu_.has_tfsf && has_tfsf_gpu_) {
				std::cout << "  TFSF debug | abs(E)=" << std::scientific << std::setprecision(3)
				          << tfsf_abs_sum_e
				          << " abs(H)=" << tfsf_abs_sum_h
				          << " max(E)=" << tfsf_max_e
				          << " max(H)=" << tfsf_max_h
				          << '\n';
			}
			std::cout << "  out  debug | abs(out)=" << out_abs_sum
			          << " max(out)=" << out_max
			          << " stale(out)=" << out_abs_sum_stale
			          << " abs(elem_out)=" << elem_out_abs_sum
			          << " max(elem_out)=" << elem_out_max
			          << '\n' << std::defaultfloat;
			std::cout << std::flush;
		}

		acc_total = acc_exchange = acc_jumps = acc_elements = 0;
		acc_promoted = acc_aux = acc_sgbc = 0;
		acc_count = 0;
		lastPrintTime = currentTime;
	}
#endif
}

void hesthaven_add_curved_gpu(HesthavenGPUData& gpu,
                              mfem::SparseMatrix& curved_merged,
                              const std::array<mfem::ParGridFunction, 3>& eOld,
                              const std::array<mfem::ParGridFunction, 3>& hOld,
                              mfem::Vector& out,
                              int ndofs,
                              int nbr_size)
{
	const int blockSize = ndofs + nbr_size;
	const int fullSize = gpu.curved_full_size;

	// Assemble ext_in on host (same path as CPU Mult) so halo data matches jumps.
	mfem::Vector ext_in(fullSize);
	ext_in.UseDevice(false);
	ext_in = 0.0;
	for (int d = maxwell::X; d <= maxwell::Z; ++d) {
		ext_in.SetVector(eOld[d],       d      * blockSize);
		ext_in.SetVector(hOld[d],  (3 + d) * blockSize);
		if (nbr_size > 0) {
			ext_in.SetVector(eOld[d].FaceNbrData(),      d      * blockSize + ndofs);
			ext_in.SetVector(hOld[d].FaceNbrData(), (3 + d) * blockSize + ndofs);
		}
	}

	std::memcpy(gpu.d_ext_in.HostWrite(), ext_in.HostRead(),
	            static_cast<size_t>(fullSize) * sizeof(double));
	(void)gpu.d_ext_in.Read();

	gpu.d_curved_y = 0.0;
	curved_merged.AddMult(gpu.d_ext_in, gpu.d_curved_y, 1.0);

	// GPU SpMV done; merge curved increments on host (proven coupling path).
	const double* y_host = gpu.d_curved_y.HostRead();
	const int* map_host = gpu.d_curved_row_to_out.HostRead();
	(void)out.HostRead();
	double* out_host = out.HostWrite();
	for (int k = 0; k < gpu.n_curved_rows; ++k) {
		out_host[map_host[k]] += y_host[k];
	}
	(void)out.Read();
}

#endif

} // namespace maxwell