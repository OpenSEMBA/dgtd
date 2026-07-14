#include "HesthavenEvolution.h"

#include "components/DGOperatorFactory.h"
#include "solver/SourcesManager.h"
#ifdef SEMBA_DGTD_ENABLE_CUDA
#include "GlobalEvolution.h"
#endif

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
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
}

} // namespace
#endif

using MatricesSet = std::set<DynamicMatrix, MatrixCompareLessThan>;

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
					// Sign matches GlobalEvolution (out -= TFSFOp * planewave).
					out.e_[d][mapBSF[m][v]] += fields[E][d];
					out.h_[d][mapBSF[m][v]] += fields[H][d];
					out.e_[d][mapBTF[m][v]] -= fields[E][d];
					out.h_[d][mapBTF[m][v]] -= fields[H][d];
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
					e_jumps[d * jumps_size + mapBSF[m][v]] += fields[E][d];
					h_jumps[d * jumps_size + mapBSF[m][v]] += fields[H][d];
					e_jumps[d * jumps_size + mapBTF[m][v]] -= fields[E][d];
					h_jumps[d * jumps_size + mapBTF[m][v]] -= fields[H][d];
				}
			}
		}
	}

}

void HesthavenEvolution::addCurvedElementContributions(Vector& out) const
{
	if (hestElemCurvedStorage_.empty()) {
		return;
	}

	const int ndofs = fes_.GetNDofs();
	const int nbrDofs = fes_.num_face_nbr_dofs;

#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (gpu_.curved_initialized && curved_merged_matrix_ &&
	    mfem::Device::Allows(mfem::Backend::CUDA) && out.UseDevice()) {
		hesthaven_add_curved_gpu(gpu_, *curved_merged_matrix_, eOld_, hOld_, out,
		                         ndofs, nbrDofs);
		return;
	}
#endif

	const int blockSize = ndofs + nbrDofs;
	const int localSize = numberOfFieldComponents * numberOfMaxDimensions * ndofs;
	const int fullSize = numberOfFieldComponents * numberOfMaxDimensions * blockSize;
	const int ncomp = numberOfFieldComponents * numberOfMaxDimensions;

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

	Vector curved_out(localSize);
	curved_out.UseDevice(false);

#ifdef SEMBA_DGTD_ENABLE_CUDA
	const bool out_is_device =
	    mfem::Device::Allows(mfem::Backend::CUDA) && out.UseDevice();
	if (out_is_device) {
		(void)out.HostRead();
	}
#endif

	double* out_data = out.HostWrite();
	for (const auto& curved : hestElemCurvedStorage_) {
		curved_out = 0.0;
		curved.matrix.AddMult(ext_in, curved_out);
		const double* co = curved_out.HostRead();
		for (int ii = 0; ii < curved.dofs.Size(); ++ii) {
			const int dof = curved.dofs[ii];
			for (int c = 0; c < ncomp; ++c) {
				out_data[c * ndofs + dof] += co[c * ndofs + dof];
			}
		}
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

void HesthavenEvolution::storeDirectionalMatrices(FiniteElementSpace& subFES, const DynamicMatrix& refInvMass, HesthavenElement& hestElem)
{
	// Single-element submesh: keep assembly serial on this rank (MPI_COMM_SELF).
	Model model(*subFES.GetMesh(), GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo{}, nullptr, MPI_COMM_SELF);
	Probes probes;
	ProblemDescription pd(model, probes, srcmngr_.sources, opts_);
	DGOperatorFactory<FiniteElementSpace> dgops(pd, subFES);
	for (int d = X; d <= Z; d++) {
		auto denseMat = dgops.buildDerivativeSubOperator<BilinearForm>(d)->SpMat().ToDenseMatrix();
		DynamicMatrix dirMat = refInvMass * toEigen(*denseMat) * getReferenceVolume(hestElem.type) / hestElem.vol;
		delete denseMat;
		StorageIterator it = matrixStorage_.find(dirMat);
		if (it == matrixStorage_.end()) {
			matrixStorage_.insert(dirMat);
			it = matrixStorage_.find(dirMat);
			hestElem.dir[d] = &(*it);
		}
		else {
			hestElem.dir[d] = &(*it);
		} 
	}
}

void storeFaceInformation(FiniteElementSpace& subFES, HesthavenElement& hestElem)
{

	int numFaces, numNodesAtFace;
	subFES.GetMesh()->SetCurvature(1);
	const auto& dim = subFES.GetMesh()->Dimension();
	dim == 2 ? numFaces = subFES.GetMesh()->GetNEdges() : numFaces = subFES.GetMesh()->GetNFaces();
	dim == 2 ? numNodesAtFace = numNodesAtFace = subFES.FEColl()->GetOrder() + 1 : numNodesAtFace = getFaceNodeNumByGeomType(subFES);
	initNormalVectors(hestElem, numFaces * numNodesAtFace);
	initFscale(hestElem, numFaces * numNodesAtFace);
	auto J{ subFES.GetMesh()->GetElementTransformation(0)->Weight() };

	ElementTransformation* faceTrans;
	for (int f = 0; f < numFaces; f++) {

		Vector normal(dim);
		dim == 2 ? faceTrans = subFES.GetMesh()->GetEdgeTransformation(f) : faceTrans = subFES.GetMesh()->GetFaceTransformation(f);
		auto ir{ &IntRules.Get(faceTrans->GetGeometryType(), faceTrans->OrderW() + 2 * faceTrans->Order()) };
		faceTrans->SetIntPoint(&ir->IntPoint(0)); // Face is always order 1, thus we don't need more than one point to calculate the normal at the face.
		CalcOrtho(faceTrans->Jacobian(), normal);
		auto sJ{ faceTrans->Weight() };

		for (auto b= 0; b < numNodesAtFace; b++) { //hesthaven requires normals to be stored once per node at face
			hestElem.normals[X][f * numNodesAtFace + b] = normal[0] / sJ;
			hestElem.fscale[f * numNodesAtFace + b] = sJ * 2.0 / J; //likewise for fscale, surface per volume ratio per node at face
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
	applyBdrBC(bdrMaps.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);

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
	applyBdrBC(bdrMaps.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);

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

HesthavenEvolution::HesthavenEvolution(ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts)
{

	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);

	const auto* cmesh = &model_.getConstMesh();
	ParMesh& mesh = model_.getMesh();
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
	auto attMap{ mapOriginalAttributes(mesh) };

	positions_ = buildDoFPositions(fes_);
	auto elemOrderList = initCurvedAndLinearElementsLists(fes_, positions_);
	linearElements_ = elemOrderList.first;
	curvedElements_ = elemOrderList.second;

	checkForTFSFInCurvedElements();

	hestElemLinearStorage_.resize(linearElements_.Size());

	bool allElementsSameGeomType = true; 
	{
		const auto firstElemGeomType = cmesh->GetElementGeometry(0);
		for (auto e= 0; e < cmesh->GetNE(); e++)
		{
			if (firstElemGeomType != cmesh->GetElementGeometry(e))
			{
				allElementsSameGeomType = false;
			}
		}
	}

	DynamicMatrix refInvMass;
	if (allElementsSameGeomType)
	{
		refInvMass = assembleHesthavenRefElemInvMassMatrix(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
		refLIFT_ = refInvMass * assembleHesthavenRefElemEmat(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
	}

	hestElemLinearStorage_.resize(linearElements_.Size());

	const int n_assemble = linearElements_.Size();
	if (Mpi::WorldRank() == 0) {
		std::cout << "Hesthaven init: assembling " << n_assemble << " linear elements...\n"
		          << std::flush;
	}

	std::atomic<int> assembled_count{0};
	const int progress_stride = std::max(1, n_assemble / 16);

	auto assembleOneElement = [&](int e, Mesh& work_mesh) {
		HesthavenElement hestElem;
		hestElem.id = linearElements_[e];
		hestElem.type = cmesh->GetElementType(linearElements_[e]);

		work_mesh.SetAttribute(linearElements_[e], hesthavenMeshingTag);
		hestElem.vol = work_mesh.GetElementVolume(linearElements_[e]);
		auto sm{ SubMesh::CreateFromDomain(work_mesh, elementMarker) };
		restoreOriginalAttributesAfterSubMeshing(linearElements_[e], work_mesh, attMap);
		FiniteElementSpace subFES(&sm, fec);

		sm.bdr_attributes.SetSize(subFES.GetNF());
		for (auto f = 0; f < subFES.GetNF(); f++) {
			sm.bdr_attributes[f] = f + 1;
			sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
		}

#ifdef SEMBA_DGTD_ENABLE_OPENMP
		#pragma omp critical(hesthaven_matrix_dedup)
#endif
		storeDirectionalMatrices(subFES, refInvMass, hestElem);
		storeFaceInformation(subFES, hestElem);
		hestElemLinearStorage_[e] = std::move(hestElem);

		const int done = ++assembled_count;
		if (Mpi::WorldRank() == 0 && (done % progress_stride == 0 || done == n_assemble)) {
			std::cout << "Hesthaven init: " << done << "/" << n_assemble << " elements\n"
			          << std::flush;
		}
	};

#ifdef SEMBA_DGTD_ENABLE_OPENMP
	bool use_parallel_assembly = (n_assemble > 32);
#ifdef SEMBA_DGTD_ENABLE_CUDA
	// MFEM mesh/DG assembly is not thread-safe with the CUDA device backend enabled.
	if (mfem::Device::Allows(mfem::Backend::CUDA)) {
		use_parallel_assembly = false;
	}
#endif
	if (use_parallel_assembly) {
		#pragma omp parallel
		{
			Mesh work_mesh;
			#pragma omp for schedule(dynamic, 64)
			for (int e = 0; e < linearElements_.Size(); e++) {
				#pragma omp critical(hesthaven_mesh_copy)
				{
					work_mesh = Mesh(mesh);
				}
				assembleOneElement(e, work_mesh);
			}
		}
	} else
#endif
	{
		for (int e = 0; e < linearElements_.Size(); e++) {
			Mesh work_mesh(mesh);
			assembleOneElement(e, work_mesh);
		}
	}

	if (curvedElements_.size()) {
		Probes probes;
		ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
		DGOperatorFactory<ParFiniteElementSpace> dgops(pd, fes_);
		auto global = dgops.buildGlobalOperator();

		const int ndofs = fes_.GetNDofs();
		const int nbrDofs = fes_.num_face_nbr_dofs;
		const int blockSize = ndofs + nbrDofs;
		const int localSize = numberOfFieldComponents * numberOfMaxDimensions * ndofs;
		const int fullSize = numberOfFieldComponents * numberOfMaxDimensions * blockSize;

		for (const auto& [e, dofs]: curvedElements_) {
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

	connectivity_.emplace(model, fes_);

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

#ifdef SEMBA_DGTD_ENABLE_CUDA
	if (linearElements_.Size() > 0 && refLIFT_.size() > 0) {
		initGPUData();
	}
	if (!hestElemCurvedStorage_.empty()) {
		initGPUCurvedData();
	}
#endif

#ifdef SHOW_TIMER_INFORMATION
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
	mfem::StopWatch timerTotal, timerExchange, timerJumps, timerElements, timerCurved;
	timerTotal.Start();
	double tfsf_abs_sum_e = 0.0;
	double tfsf_abs_sum_h = 0.0;
	double tfsf_max_e = 0.0;
	double tfsf_max_h = 0.0;
	double out_abs_sum = 0.0;
	double out_max = 0.0;
#endif

	double alpha = opts_.alpha;
	(void)in.HostRead();

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

	auto scatterElem = [&](const Eigen::VectorXd& data, ElementId id,
	                       Eigen::Map<Eigen::VectorXd>& out_map) {
		Array<int> dofs;
		fes_.GetElementDofs(id, dofs);
		for (int i = 0; i < dofs.Size(); ++i) {
			out_map[dofs[i]] += data[i];
		}
	};

	auto jumps{ HesthavenFields(connectivity_->global.size()) };
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

#ifdef SEMBA_DGTD_ENABLE_OPENMP
	#pragma omp parallel for schedule(static) if(linearElements_.Size() > 32)
#endif
	for (int e = 0 ; e < linearElements_.Size(); e++) {
		auto elemFluxSize{ hestElemLinearStorage_[e].fscale.size() };

		const auto jumpsElem{ HesthavenElementJumps(jumps, hestElemLinearStorage_[e].id, elemFluxSize) };
		const auto fieldsElem{ FieldsElementMaps(in, fes_, hestElemLinearStorage_[e].id) };

		Eigen::VectorXd ndotdH(jumpsElem.h_[X].size()), ndotdE(jumpsElem.e_[X].size());
		ndotdH.setZero(); ndotdE.setZero();

		for (int d = X; d <= Z; d++) {
			ndotdH += hestElemLinearStorage_[e].normals[d].asDiagonal() * jumpsElem.h_[d];
			ndotdE += hestElemLinearStorage_[e].normals[d].asDiagonal() * jumpsElem.e_[d];
		}

		HesthavenFields elemFlux(elemFluxSize);

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const Eigen::VectorXd& norx = hestElemLinearStorage_[e].normals[x];
			const Eigen::VectorXd& nory = hestElemLinearStorage_[e].normals[y];
			const Eigen::VectorXd& norz = hestElemLinearStorage_[e].normals[z];

			elemFlux.h_[x] = -1.0 * nory.asDiagonal() * jumpsElem.e_[z] +       norz.asDiagonal() * jumpsElem.e_[y] + alpha * (jumpsElem.h_[x] - ndotdH.asDiagonal() * norx);
			elemFlux.e_[x] =        nory.asDiagonal() * jumpsElem.h_[z] - 1.0 * norz.asDiagonal() * jumpsElem.h_[y] + alpha * (jumpsElem.e_[x] - ndotdE.asDiagonal() * norx);

		}

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const DynamicMatrix& dir1 = *hestElemLinearStorage_[e].dir[y];
			const DynamicMatrix& dir2 = *hestElemLinearStorage_[e].dir[z];

			const Eigen::VectorXd& hResult = -1.0 * dir1 * fieldsElem.e_[z] +        dir2 * fieldsElem.e_[y] + applyLIFT(hestElemLinearStorage_[e].fscale, elemFlux.h_[x]);
			const Eigen::VectorXd& eResult =        dir1 * fieldsElem.h_[z] - 1.0 *  dir2 * fieldsElem.h_[y] + applyLIFT(hestElemLinearStorage_[e].fscale, elemFlux.e_[x]);

			scatterElem(hResult, hestElemLinearStorage_[e].id, outMaps.h_[x]);
			scatterElem(eResult, hestElemLinearStorage_[e].id, outMaps.e_[x]);
		}

	}
#ifdef SHOW_TIMER_INFORMATION
	timerElements.Stop();
	timerCurved.Start();
#endif

	addCurvedElementContributions(out);

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
	timerCurved.Stop();
	timerTotal.Stop();

	static double acc_total = 0, acc_exchange = 0, acc_jumps = 0, acc_elements = 0, acc_curved = 0;
	static int acc_count = 0;
	static double lastPrintTime = -1.0;
	constexpr double printInterval = 0.05;
	constexpr int printCadence = 50;

	acc_total    += timerTotal.RealTime()    * 1000.0;
	acc_exchange += timerExchange.RealTime() * 1000.0;
	acc_jumps    += timerJumps.RealTime()    * 1000.0;
	acc_elements += timerElements.RealTime() * 1000.0;
	acc_curved   += timerCurved.RealTime()   * 1000.0;
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
		double local_avg[5] = {
			acc_total    / acc_count,
			acc_exchange / acc_count,
			acc_jumps    / acc_count,
			acc_elements / acc_count,
			acc_curved   / acc_count
		};

		std::vector<double> all_avg(5 * P);
		if (P > 1) {
			MPI_Gather(local_avg, 5, MPI_DOUBLE,
			           all_avg.data(), 5, MPI_DOUBLE,
			           0, fes_.GetParMesh()->GetComm());
		} else {
			std::copy(local_avg, local_avg + 5, all_avg.data());
		}

		if (Mpi::WorldRank() == 0) {
			std::cout << "\n[Hesthaven Mult timing] t=" << std::fixed << std::setprecision(4)
			          << currentTime << "  (avg of " << acc_count << " calls, ms/call)\n"
			          << std::defaultfloat;
			std::cout << "  Rank  | total   exchange    jumps  elements   curved  | bottleneck\n"
			          << "  ------+--------------------------------------------------+-----------\n";
			for (int r = 0; r < P; ++r) {
				const double t_total    = all_avg[r * 5 + 0];
				const double t_exchange = all_avg[r * 5 + 1];
				const double t_jumps    = all_avg[r * 5 + 2];
				const double t_elements = all_avg[r * 5 + 3];
				const double t_curved   = all_avg[r * 5 + 4];
				const char* bottleneck = "elements";
				if (t_exchange > 0.7 * t_total) {
					bottleneck = "MPI wait";
				} else if (t_jumps > 0.5 * t_total) {
					bottleneck = "jumps/BC";
				} else if (t_curved > 0.3 * t_total) {
					bottleneck = "curved";
				}

				char buf[256];
				std::snprintf(buf, sizeof(buf),
				              "  %4d  | %6.2f  %8.2f  %7.2f  %8.2f  %7.3f  | %s",
				              r, t_total, t_exchange, t_jumps, t_elements, t_curved, bottleneck);
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

		acc_total = acc_exchange = acc_jumps = acc_elements = acc_curved = 0;
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
		Array<int> el_dofs;
		fes_.GetElementDofs(he.id, el_dofs);
		for (int i = 0; i < el_dofs.Size(); ++i) {
			elem_dofs_host[e * gpu_.ndof_el + i] = el_dofs[i];
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

	hesthaven_sync_gpu_static_data(gpu_);

	gpu_.initialized = true;
}

void HesthavenEvolution::initGPUCurvedData()
{
	if (!mfem::Device::Allows(mfem::Backend::CUDA)) {
		return;
	}
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
	gpu_.d_curved_row_to_out.SetSize(n_rows);
	gpu_.d_curved_row_to_out.GetMemory().UseDevice(true);

	int row = 0;
	int* row_to_out = gpu_.d_curved_row_to_out.HostWrite();
	for (const auto& curved : hestElemCurvedStorage_) {
		for (int ii = 0; ii < curved.dofs.Size(); ++ii) {
			const int dof = curved.dofs[ii];
			for (int c = 0; c < ncomp; ++c) {
				const int global_row = c * ndofs + dof;
				Array<int> cols;
				Vector vals;
				curved.matrix.GetRow(global_row, cols, vals);
				curved_merged_matrix_->SetRow(row, cols, vals);
				row_to_out[row] = global_row;
				++row;
			}
		}
	}
	curved_merged_matrix_->Finalize();

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
	appendTrueBC(bdr.SMA, -1.0 / opts_.alpha, -1.0 / opts_.alpha);

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
	mfem::StopWatch timerTotal, timerExchange, timerJumps, timerElements, timerCurved;
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

	out.SetSize(Height());
	out.UseDevice(true);
	out = 0.0;

#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Start();
#endif
	exchangeFieldData(in, false);
#ifdef SHOW_TIMER_INFORMATION
	timerExchange.Stop();
	timerJumps.Start();
#endif

	// Jumps, BC, and TFSF on host; MPI face-neighbor data arrives in host buffers.
	{
		(void)in.HostRead();
		auto jumps{ HesthavenFields(gpu_.jumps_size) };
		computeJumps(in, jumps);
		const FieldsInputMaps fieldsIn(in, fes_);
		applyBoundaryConditionsToNodes(connectivity_->boundary, fieldsIn, jumps);
		evaluateTFSF(jumps);

		const int js = gpu_.jumps_size;
		double* jumps_e_host = gpu_.d_jumps_e.HostWrite();
		double* jumps_h_host = gpu_.d_jumps_h.HostWrite();
		for (int d = X; d <= Z; ++d) {
			std::memcpy(jumps_e_host + d * js,
			            jumps.e_[d].data(),
			            static_cast<size_t>(js) * sizeof(double));
			std::memcpy(jumps_h_host + d * js,
			            jumps.h_[d].data(),
			            static_cast<size_t>(js) * sizeof(double));
		}
		// HostWrite above; push to device only (ReadWrite would pull stale device data).
		(void)gpu_.d_jumps_e.Read();
		(void)gpu_.d_jumps_h.Read();

#ifdef SEMBA_DGTD_ENABLE_CUDA
		// Upload MPI halos to device only after host-side jump assembly is done.
		if (mfem::Device::Allows(mfem::Backend::CUDA) && fes_.num_face_nbr_dofs > 0) {
			sync_cuda_face_nbr_halos(eOld_, hOld_);
		}
#endif

#ifdef SHOW_TIMER_INFORMATION
		if (gpu_.has_tfsf) {
			for (int d = X; d <= Z; ++d) {
				for (int i = 0; i < js; ++i) {
					const double ae = std::abs(jumps.e_[d][i]);
					const double ah = std::abs(jumps.h_[d][i]);
					tfsf_abs_sum_e += ae;
					tfsf_abs_sum_h += ah;
					tfsf_max_e = std::max(tfsf_max_e, ae);
					tfsf_max_h = std::max(tfsf_max_h, ah);
				}
			}
		}
#endif
	}

#ifdef SHOW_TIMER_INFORMATION
	timerJumps.Stop();
	timerElements.Start();
#endif

	if (linearElements_.Size() > 0) {
		// Jump path used HostRead(); push the same state to device for the element kernel.
		const_cast<mfem::Vector&>(in).Read();
		hesthaven_mult_gpu(gpu_, opts_.alpha, in, eOld_, hOld_, out, ndofs);
	}

#ifdef SHOW_TIMER_INFORMATION
	timerElements.Stop();
	timerCurved.Start();
#endif

	addCurvedElementContributions(out);

#ifdef SHOW_TIMER_INFORMATION
	{
		const double* stale_out = out.HostRead();
		for (int i = 0; i < out.Size(); ++i) {
			const double a = std::abs(stale_out[i]);
			out_abs_sum_stale += a;
			out_max_stale = std::max(out_max_stale, a);
		}
	}
#endif

	// Device-authoritative 'out' must be synced before host-side debug/export reads.
	(void)out.Read();

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
	timerCurved.Stop();
	timerTotal.Stop();

	static double acc_total = 0, acc_exchange = 0, acc_jumps = 0, acc_elements = 0, acc_curved = 0;
	static int acc_count = 0;
	static double lastPrintTime = -1.0;
	constexpr double printInterval = 0.05;
	constexpr int printCadence = 50;

	acc_total    += timerTotal.RealTime()    * 1000.0;
	acc_exchange += timerExchange.RealTime() * 1000.0;
	acc_jumps    += timerJumps.RealTime()    * 1000.0;
	acc_elements += timerElements.RealTime() * 1000.0;
	acc_curved   += timerCurved.RealTime()   * 1000.0;
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
		double local_avg[5] = {
			acc_total    / acc_count,
			acc_exchange / acc_count,
			acc_jumps    / acc_count,
			acc_elements / acc_count,
			acc_curved   / acc_count
		};

		std::vector<double> all_avg(5 * P);
		if (P > 1) {
			MPI_Gather(local_avg, 5, MPI_DOUBLE,
			           all_avg.data(), 5, MPI_DOUBLE,
			           0, fes_.GetParMesh()->GetComm());
		} else {
			std::copy(local_avg, local_avg + 5, all_avg.data());
		}

		if (Mpi::WorldRank() == 0) {
			std::cout << "\n[Hesthaven Mult timing] t=" << std::fixed << std::setprecision(4)
			          << currentTime << "  (avg of " << acc_count << " calls, ms/call, GPU)\n"
			          << std::defaultfloat;
			std::cout << "  Rank  | total   exchange    jumps  elements   curved  | bottleneck\n"
			          << "  ------+--------------------------------------------------+-----------\n";
			for (int r = 0; r < P; ++r) {
				const double t_total    = all_avg[r * 5 + 0];
				const double t_exchange = all_avg[r * 5 + 1];
				const double t_jumps    = all_avg[r * 5 + 2];
				const double t_elements = all_avg[r * 5 + 3];
				const double t_curved   = all_avg[r * 5 + 4];
				const char* bottleneck = "elements";
				if (t_exchange > 0.7 * t_total) {
					bottleneck = "MPI wait";
				} else if (t_jumps > 0.5 * t_total) {
					bottleneck = "jumps/BC";
				} else if (t_curved > 0.3 * t_total) {
					bottleneck = "curved";
				}

				char buf[256];
				std::snprintf(buf, sizeof(buf),
				              "  %4d  | %6.2f  %8.2f  %7.2f  %8.2f  %7.3f  | %s",
				              r, t_total, t_exchange, t_jumps, t_elements, t_curved, bottleneck);
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

		acc_total = acc_exchange = acc_jumps = acc_elements = acc_curved = 0;
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