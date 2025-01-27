#include "Evolution.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

void evalConductivity(const Vector& cond, const Vector& in, Vector& out)
{
	for (auto v{ 0 }; v < cond.Size(); v++) {
		out[v] -= cond[v] * in[v];
	}
}

void changeSignOfFieldGridFuncs(FieldGridFuncs& gfs)
{
	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			gfs[f][d] *= -1.0;
		}
	}
}

const FieldGridFuncs evalTimeVarFunction(const Time time, SourcesManager& sm)
{
	auto res{ sm.evalTimeVarField(time, sm.getGlobalTFSFSpace()) };
	auto func_g_sf = res;
	sm.markDoFSforTForSF(res, true);
	{
		if (sm.getTFSFSubMesher().getSFSubMesh() != NULL) {
			sm.markDoFSforTForSF(func_g_sf, false);
			for (int f : {E, H}) {
				for (int x{ 0 }; x <= Z; x++) {
					res[f][x] -= func_g_sf[f][x];
					res[f][x] *= 0.5;
				}
			}
		}
	}
	return res;
}

Evolution::Evolution(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
#ifdef SHOW_TIMER_INFORMATION
	auto startTime{ std::chrono::high_resolution_clock::now() };
#endif


#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "---------OPERATOR ASSEMBLY INFORMATION----------" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
#endif

	if (model_.getTotalFieldScatteredFieldToMarker().find(BdrCond::TotalFieldIn) != model_.getTotalFieldScatteredFieldToMarker().end()) {
		srcmngr_.initTFSFPreReqs(model_.getConstMesh(), model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn));
		auto globalTFSFfes{ srcmngr_.getGlobalTFSFSpace() };
		Model modelGlobal = Model(*globalTFSFfes->GetMesh(), GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Assembling TFSF Inverse Mass Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MInvTFSF_[f] = buildInverseMassMatrix(f, modelGlobal, *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono:: >
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MP_GTFSF_[f] = buildByMult(*MInvTFSF_[f], *buildZeroNormalOperator(f, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					MFN_GTFSF_[f][f2][d] = buildByMult(*MInvTFSF_[f], *buildOneNormalOperator(f2, {d}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
				}
			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						MFNN_GTFSF_[f][f2][d][d2] = buildByMult(*MInvTFSF_[f], *buildTwoNormalOperator(f2, {d, d2}, modelGlobal, *globalTFSFfes, opts_), *globalTFSFfes);
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		MInv_[f] = buildInverseMassMatrix(f, model_, fes_);
	}

	if (model_.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			MPB_[f] = buildByMult(*MInv_[f], *buildZeroNormalIBFIOperator(f, model_, fes_, opts_), fes_);
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					MFNB_[f][f2][d] = buildByMult(*MInv_[f], *buildOneNormalIBFIOperator(f2, { d }, model_, fes_, opts_), fes_);
				} 
			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
#endif

		for (auto f : { E, H }) {
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						MFNNB_[f][f2][d][d2] = buildByMult(*MInv_[f], *buildTwoNormalIBFIOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
					}
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			MS_[f][d] = buildByMult(*MInv_[f], *buildDerivativeOperator(d, fes_), fes_);
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		MP_[f] = buildByMult(*MInv_[f], *buildZeroNormalOperator(f, model_, fes_, opts_), fes_);
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			for (auto f2 : { E, H }) {
				MFN_[f][f2][d] = buildByMult(*MInv_[f], *buildOneNormalOperator(f2, { d }, model_, fes_, opts_), fes_);
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
#endif

	for (auto f : { E, H }) {
		for (auto d{ X }; d <= Z; d++) {
			for (auto f2 : { E, H }) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					MFNN_[f][f2][d][d2] = buildByMult(*MInv_[f], *buildTwoNormalOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
				}
			}
		}
	}

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "Elapsed time (ms): " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
	std::cout << "Operator assembly finished" << std::endl;
	std::cout << std::endl;
#endif

	CND_ = buildConductivityCoefficients(model_, fes_);
 }

void Evolution::Mult(const Vector& in, Vector& out) const
{
	const auto& dim{ fes_.GetMesh()->Dimension() };

	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		eNew[d].SetSpace(&fes_);
		hNew[d].SetSpace(&fes_);
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;
		
		evalConductivity(CND_, eOld[x], eNew[x]);

		//Centered
		MS_[H][y]		 ->AddMult(eOld[z], hNew[x],-1.0);
		MS_[H][z]		 ->AddMult(eOld[y], hNew[x]);
		MS_[E][y]		 ->AddMult(hOld[z], eNew[x]);
		MS_[E][z]		 ->AddMult(hOld[y], eNew[x],-1.0);
		
		MFN_[H][E][y]  	 ->AddMult(eOld[z], hNew[x], 1.0);
		MFN_[H][E][z]	 ->AddMult(eOld[y], hNew[x],-1.0);
		MFN_[E][H][y]  	 ->AddMult(hOld[z], eNew[x],-1.0);
		MFN_[E][H][z]	 ->AddMult(hOld[y], eNew[x], 1.0);

		if (opts_.fluxType == FluxType::Upwind) {

			MFNN_[H][H][X][x]->AddMult(hOld[X], hNew[x], 1.0);
			MFNN_[H][H][Y][x]->AddMult(hOld[Y], hNew[x], 1.0);
			MFNN_[H][H][Z][x]->AddMult(hOld[Z], hNew[x], 1.0);
			MP_[H]			 ->AddMult(hOld[x], hNew[x],-1.0);
			
			MFNN_[E][E][Y][x]->AddMult(eOld[Y], eNew[x], 1.0);
			MFNN_[E][E][X][x]->AddMult(eOld[X], eNew[x], 1.0);
			MFNN_[E][E][Z][x]->AddMult(eOld[Z], eNew[x], 1.0);
			MP_[E]			 ->AddMult(eOld[x], eNew[x],-1.0);
		}

		if (model_.getInteriorBoundaryToMarker().size() != 0) {
			
			MFNB_[H][E][y]->AddMult(eOld[z], hNew[x]);
			MFNB_[H][E][z]->AddMult(eOld[y], hNew[x], -1.0);
			MFNB_[E][H][y]->AddMult(hOld[z], eNew[x], -1.0);
			MFNB_[E][H][z]->AddMult(hOld[y], eNew[x]);

			if (opts_.fluxType == FluxType::Upwind) {

				MFNNB_[H][H][X][x]->AddMult(hOld[X], hNew[x]);
				MFNNB_[H][H][Y][x]->AddMult(hOld[Y], hNew[x]);
				MFNNB_[H][H][Z][x]->AddMult(hOld[Z], hNew[x]);
				MPB_[H]->AddMult(hOld[x], hNew[x], -1.0);

				MFNNB_[E][E][X][x]->AddMult(eOld[X], eNew[x]);
				MFNNB_[E][E][Y][x]->AddMult(eOld[Y], eNew[x]);
				MFNNB_[E][E][Z][x]->AddMult(eOld[Z], eNew[x]);
				MPB_[E]->AddMult(eOld[x], eNew[x], -1.0);
			}
		}
	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			
			auto func { evalTimeVarFunction(GetTime(),srcmngr_) };

			std::array<GridFunction, 3> eTemp, hTemp;

			for (int d = X; d <= Z; d++) {
				eTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				hTemp[d].SetSpace(srcmngr_.getGlobalTFSFSpace());
				eTemp[d] = 0.0;
				hTemp[d] = 0.0;
			}

			for (int x = X; x <= Z; x++) {
				int y = (x + 1) % 3;
				int z = (x + 2) % 3;

				MaxwellTransferMap eMap(eTemp[x], eNew[x]);
				MaxwellTransferMap hMap(hTemp[x], hNew[x]);

				//Centered

				MFN_GTFSF_[H][E][y]->Mult(func[E][z], hTemp[x]);
				eMap.TransferSub(hTemp[x], hNew[x]);
				MFN_GTFSF_[H][E][z]->Mult(func[E][y], hTemp[x]);
				eMap.TransferAdd(hTemp[x], hNew[x]);
				MFN_GTFSF_[E][H][y]->Mult(func[H][z], eTemp[x]);
				eMap.TransferAdd(eTemp[x], eNew[x]);
				MFN_GTFSF_[E][H][z]->Mult(func[H][y], eTemp[x]);
				eMap.TransferSub(eTemp[x], eNew[x]);

				if (opts_.fluxType == FluxType::Upwind) {
					MFNN_GTFSF_[H][H][X][x]->Mult(func[H][X], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Y][x]->Mult(func[H][Y], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MFNN_GTFSF_[H][H][Z][x]->Mult(func[H][Z], hTemp[x]);
					hMap.TransferSub(hTemp[x], hNew[x]);
					MP_GTFSF_[H]           ->Mult(func[H][x], hTemp[x]);
					hMap.TransferAdd(hTemp[x], hNew[x]);

					MFNN_GTFSF_[E][E][X][x]->Mult(func[E][X], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Y][x]->Mult(func[E][Y], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MFNN_GTFSF_[E][E][Z][x]->Mult(func[E][Z], eTemp[x]);
					eMap.TransferSub(eTemp[x], eNew[x]);
					MP_GTFSF_[E]           ->Mult(func[E][x], eTemp[x]);
					eMap.TransferAdd(eTemp[x], eNew[x]);
				}

			}
		}
	}
}

DynamicMatrix assembleInverseMassMatrix(FiniteElementSpace& fes)
{
	BilinearForm bf(&fes);
	ConstantCoefficient one(1.0);
	bf.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	bf.Assemble();
	bf.Finalize();

	return toEigen(*bf.SpMat().ToDenseMatrix());
}

DynamicMatrix getReferenceInverseMassMatrix(const Element::Type elType, const int order, const int dimension)
{
	std::unique_ptr<Mesh> m;
	switch (dimension) {
	case 2:
		m = std::make_unique<Mesh>(Mesh::MakeCartesian2D(1, 1, elType));
		break;
	case 3:
		m = std::make_unique<Mesh>(Mesh::MakeCartesian3D(1, 1, 1, elType));
		break;
	default:
		throw std::runtime_error("Hesthaven Evolution Operator only supports dimensions 2 or 3.");
	}
	auto fec{ L2_FECollection(order, dimension, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(m.get(), &fec)};
	auto mass_mat{ assembleInverseMassMatrix(fes) };

	DynamicMatrix res = getElementMassMatrixFromGlobal(0, mass_mat);
	return res;
}

GlobalBoundary assembleGlobalBoundaryMap(const GlobalConnectivity& vmaps, BoundaryToMarker& markers, FiniteElementSpace& fes)
{
	GlobalBoundary res;

	for (auto& [bdr_cond, marker] : markers) {
		auto bf{ BilinearForm(&fes) };
		bf.AddBdrFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0), markers[bdr_cond]);
		bf.Assemble();
		bf.Finalize();
		auto bdr_matrix{ toEigen(*bf.SpMat().ToDenseMatrix()) };
		std::cout << bdr_matrix << std::endl;
		std::vector<int> nodes;
		for (auto r{ 0 }; r < bdr_matrix.rows(); r++) {
			if (bdr_matrix(r, r) != 0.0) { //These conditions would be those of a node on itself, thus we only need to check if the 'self-value' is not zero.
				auto id = std::distance(vmaps.begin(), std::find_if(vmaps.begin(), vmaps.end(), [&](const auto& pair) { return pair.first == r && pair.second == r; }));
				nodes.push_back(id);
			}
		}
		res.push_back(std::make_pair(bdr_cond, nodes));
	}

	return res;
}

GlobalInteriorBoundary assembleGlobalInteriorBoundaryMap(InteriorBoundaryToMarker& markers, FiniteElementSpace& fes)
{
	GlobalInteriorBoundary res;

	for (auto& [bdr_cond, marker] : markers) {
		auto bf{ BilinearForm(&fes) };
		bf.AddInternalBoundaryFaceIntegrator(new mfemExtension::HesthavenFluxIntegrator(1.0), markers[bdr_cond]);
		bf.Assemble();
		bf.Finalize();
		auto int_bdr_matrix{ toEigen(*bf.SpMat().ToDenseMatrix()) };
		std::vector<int> nodes;
		for (auto r{ 0 }; r < int_bdr_matrix.rows(); r++) {
			if (int_bdr_matrix(r, r) != 0.0) {
				nodes.push_back(r);
			}
		}
		res.push_back(std::make_pair(bdr_cond, nodes));
	}

	return res;
}

const Eigen::VectorXd applyScalingFactors(const HesthavenElement& hestElem, const Eigen::VectorXd& flux)
{
	Eigen::VectorXd res(hestElem.invmass->rows());
	const auto numFaces{ getNumFaces(hestElem.geom) };
	const auto& invmass = *hestElem.invmass;
	const auto& fscale = hestElem.fscale.asDiagonal();
	const auto ematRows = hestElem.emat[0]->rows();
	const auto ematCols = hestElem.emat[0]->cols();
	DynamicMatrix emat(ematRows, ematCols * numFaces);
	for (auto f{ 0 }; f < numFaces; f++) {
		emat.block(0, f * ematCols, ematRows, ematCols) = hestElem.emat[f]->block(0, 0, ematRows, ematCols);
	}
	res = invmass * emat * fscale * flux / 2.0;
	return res;
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

void HesthavenEvolution::assembleTFSFConnectivity(const DynamicMatrix& matrix, FaceElementTransformations* faceTrans, double faceOri)
{
	for (auto r{ 0 }; r < matrix.rows() / 2; r++) {
		for (auto c{ 0 }; c < matrix.cols(); c++) {
			if (matrix(r, c) != 0.0 && c != matrix.cols() - 1) {
				for (auto c2{ c + 1 }; c2 < matrix.cols(); c2++) {
					if (matrix(r, c2) != 0.0) {
						auto elem1GlobalNodeID = (faceTrans->Elem1No * int(matrix.cols()) / 2) + c;
						auto elem2GlobalNodeID = (faceTrans->Elem2No * int(matrix.cols()) / 2) + c2 - (int(matrix.cols()) / 2);
						TFSFSigns signs;
						faceOri >= 0.0 ? signs = std::make_pair(-1.0, 1.0) : signs = std::make_pair(-1.0, 1.0); // Face Ori >= 0.0 means elem1 is SF, elem2 is TF.
						tfsf_connectivity_.push_back(std::make_pair(std::make_pair(elem1GlobalNodeID, elem2GlobalNodeID), signs));
					}
				}
			}
		}
	}
}

void HesthavenEvolution::evaluateTFSF(HesthavenFields& jumps) const
{
	auto func{ evalTimeVarFunction(GetTime(),srcmngr_) };
	for (auto i{ 0 }; i < tfsf_connectivity_.size(); i++) {
		auto it = std::find(connectivity_.begin(), connectivity_.end(), tfsf_connectivity_[i].first);
		auto connId = std::distance(std::begin(connectivity_), it);
		for (int d = X; d <= Z; d++) {
			jumps.e_[d][connId] += func[E][d][tfsf_connectivity_[i].first.first] * tfsf_connectivity_[i].second.first;
			jumps.h_[d][connId] += func[H][d][tfsf_connectivity_[i].first.second] * tfsf_connectivity_[i].second.second;
		}
	}
}

std::vector<Array<int>> assembleBdrMarkersForAllBdrElements(Mesh& bdrMesh, const int numBdrElems)
{
	std::vector<Array<int>> res(numBdrElems);
	bdrMesh.bdr_attributes.SetSize(numBdrElems);
	for (auto b{ 0 }; b < numBdrElems; b++) {
		bdrMesh.bdr_attributes[b] = b + 1;
		bdrMesh.SetBdrAttribute(b, bdrMesh.bdr_attributes[b]);
		Array<int> bdr_marker;
		bdr_marker.SetSize(numBdrElems);
		bdr_marker = 0;
		bdr_marker[b] = 1;
		res[b] = bdr_marker;
	}
	return res;
}

std::vector<NodeId> findNodesPerBdrFace(const DynamicMatrix& bdrNodeMat)
{
	std::vector<NodeId> res;
	for (auto r{ 0 }; r < bdrNodeMat.rows(); r++) {
		if (abs(bdrNodeMat(r, r)) > 1e-5) {
			for (auto c{ r }; c < bdrNodeMat.cols(); c++) {
				if (abs(bdrNodeMat(r, c)) > 1e-5) {
					res.push_back(c);
				}
			}
			break;
		}
	}
	return res;
}

std::vector<std::vector<NodeId>> assembleNodeVectorPerBdrFace(std::vector<Array<int>>& bdrNodeMarkers, Mesh& bdrMesh, const L2_FECollection* fec)
{
	auto bdrNodeFES = FiniteElementSpace(&bdrMesh, fec);
	std::vector<std::vector<NodeId>> res(bdrNodeFES.GetNBE());

	for (auto b{ 0 }; b < bdrNodeMarkers.size(); b++) {
		auto bdrNodeFinderOperator{ assembleFaceMassBilinearForm(bdrNodeFES, bdrNodeMarkers[b]) };
		auto bdrNodeMat{ toEigen(*bdrNodeFinderOperator->SpMat().ToDenseMatrix()) };
		res[b] =  findNodesPerBdrFace(bdrNodeMat);
	}

	return res;
}

std::vector<NodeId> initMapB(const GlobalConnectivity& connectivity)
{
	std::vector<NodeId> res;
	for (auto i{ 0 }; i < connectivity.size(); i++) {
		if (connectivity[i].first == connectivity[i].second) {
			res.push_back(i);
		}
	}
	return res;
}


std::vector<NodeId> initVMapB(const GlobalConnectivity& connectivity, const std::vector<NodeId>& mapB)
{
	std::vector<NodeId> res(mapB.size());
	for (auto i{ 0 }; i < mapB.size(); i++) {
		res[i] = connectivity[mapB[i]].first;
	}
	return res;
}

void HesthavenEvolution::initBdrConnectivityMaps(const std::vector<std::vector<NodeId>>& bdr2nodes)
{
	auto mapB{ initMapB(connectivity_) };
	auto vmapB{ initVMapB(connectivity_, mapB) };

	auto modelBdrMarkers = model_.getBoundaryToMarker(); //Only for true boundary types (not interior)...
	for (auto b{ 0 }; b < bdr2nodes.size(); b++) { //For each one of our bdr elements...
		for (const auto& marker : modelBdrMarkers) { //Fetch each one of our bdrMarkers in the model...
			if (marker.second[fes_.GetBdrAttribute(b) - 1] == 1) { //And check if that BdrCond is active for the specified bdrAtt of that bdr element, which means the bdr element is of that type...
				for (auto n{ 0 }; n < vmapB.size(); n++) { //Then sweep vmapB...
					if (bdr2nodes[b][0] == connectivity_[mapB[n]].first) { //To find the nodes in vmapB that are equal to the first node in our bdr element (as they are sorted when built), because it can happen...
						std::vector<NodeId> mapBToStore(bdr2nodes[b].size()); //... that we have multiple instances of a Hesthaven type node appearing twice on vmapB (i.e. corner node)...
						std::vector<NodeId> vmapBToStore(bdr2nodes[b].size()); //... so we need to save the mapB 2 vmapB pairs so the nodes are properly linked when defining boundary conditions...
						for (auto m{ 0 }; m < bdr2nodes[b].size(); m++) {
							mapBToStore[m] = mapB[n + m]; //Then make a temporary vector that we'll fill the next nodesize worth of nodes for mapB and vmapB, starting at the compared and equal node...
							vmapBToStore[m] = connectivity_[mapB[n + m]].first;
						}
						if (bdr2nodes[b] == vmapBToStore) { //As bdr2nodes[b] is already sorted when built, potVec also has to be sorted to match...
							bdr_connectivity_.mapB.push_back(mapBToStore);
							bdr_connectivity_.vmapB.push_back(std::make_pair(marker.first, vmapBToStore));  //If that's true, these nodes are of the checked marker bdrcond type, which can be PEC, PMC, etc...
							break; // And as we only support one bdrcond per face/edge/bdrwhatever, we're done for this bdr element...
						}     // As this method is completely non-dependent on dimensions or geometries, could be assumed is as generic as it could get.
					}
				}
			}
		}
	}
}

HesthavenEvolution::HesthavenEvolution(FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts)
{
	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);

	auto mesh{ Mesh(model_.getMesh()) };
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };

	connectivity_ =  assembleGlobalConnectivityMap(mesh, fec);

	auto bdrNodeMesh{ Mesh(model_.getMesh()) };
	std::vector<Array<int>> bdrNodeMarkers{ assembleBdrMarkersForAllBdrElements(bdrNodeMesh, fes_.GetNBE()) };
	auto bdr2nodes{ assembleNodeVectorPerBdrFace(bdrNodeMarkers, bdrNodeMesh, fec) };
	initBdrConnectivityMaps(bdr2nodes);
	
	//if (model.getInteriorBoundaryToMarker().size() != 0) {
	//	int_bdr_connectivity_ = assembleGlobalInteriorBoundaryMap(model.getInteriorBoundaryToMarker(), fes_);
	//}

	const auto* cmesh = &model_.getConstMesh();
	auto attMap{ mapOriginalAttributes(model_.getMesh()) };

	if (model_.getTotalFieldScatteredFieldToMarker().size() != 0) {
		auto mesh_tfsf{ Mesh(model_.getMesh()) };
		for (auto b{ 0 }; b < mesh_tfsf.GetNBE(); b++) {
			if (mesh_tfsf.GetBdrAttribute(b) == model_.getTotalFieldScatteredFieldToMarker().at(BdrCond::TotalFieldIn).Find(1) + 1) {
				auto faceOri{ calculateFaceOrientation(mesh, b) };
				auto faceTrans{ mesh_tfsf.GetInternalBdrFaceTransformations(b) };
				auto sm{ assembleInteriorFaceSubMesh(mesh, *faceTrans, attMap) };
				FiniteElementSpace smFES(&sm, fec);
				auto matrix{ assembleInteriorFluxMatrix(smFES) };
				assembleTFSFConnectivity(matrix, faceTrans, faceOri);
			}
		}
	}

	mesh = Mesh(model_.getMesh()); 

	hestElemStorage_.resize(cmesh->GetNE());

	for (auto e{ 0 }; e < cmesh->GetNE(); e++)
	{
		HesthavenElement hestElem;
		hestElem.id = e;
		hestElem.geom = cmesh->GetElementBaseGeometry(e);

		mesh.SetAttribute(e, hesthavenMeshingTag);
		auto sm = SubMesh::CreateFromDomain(mesh, elementMarker);
		restoreOriginalAttributesAfterSubMeshing(e, mesh, attMap);
		FiniteElementSpace subFES(&sm, fec);

		sm.bdr_attributes.SetSize(subFES.GetNF());
		for (auto f{ 0 }; f < subFES.GetNF(); f++) {
			sm.bdr_attributes[f] = f + 1;
			sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
		}

		auto inverseMassMatrix = assembleInverseMassMatrix(subFES);
		//auto inverseMassMatrix{ getReferenceInverseMassMatrix(sm.GetElementType(0), fec->GetOrder(), sm.Dimension())};
		StorageIterator it = matrixStorage_.find(inverseMassMatrix);
		if (it == matrixStorage_.end()) {
			matrixStorage_.insert(inverseMassMatrix);
			StorageIterator it = matrixStorage_.find(inverseMassMatrix);
			hestElem.invmass = &(*it);
		}
		else {
			hestElem.invmass = &(*it);
		}

		int numFaces;
		sm.Dimension() == 2 ? numFaces = sm.GetNEdges() : numFaces = sm.GetNFaces();

		auto boundaryMarkers = assembleBoundaryMarkers(subFES);
		for (auto f{ 0 }; f < numFaces; f++){
			auto surfaceMatrix{ assembleConnectivityFaceMassMatrix(subFES, boundaryMarkers[f]) };
			StorageIterator it = matrixStorage_.find(surfaceMatrix);
			if (it == matrixStorage_.end()) {
				matrixStorage_.insert(surfaceMatrix);
				StorageIterator it = matrixStorage_.find(surfaceMatrix);
				hestElem.emat.push_back(&(*it));
			}
			else {
				hestElem.emat.push_back(&(*it));
			}
		}

		for (auto d{ X }; d <= Z; d++) {
			auto derivativeMatrix{ toEigen(*buildDerivativeOperator(d, subFES)->SpMat().ToDenseMatrix())};
			StorageIterator it = matrixStorage_.find(derivativeMatrix);
			if (it == matrixStorage_.end()) {
				matrixStorage_.insert(derivativeMatrix);
				StorageIterator it = matrixStorage_.find(derivativeMatrix);
				hestElem.dir[d] = &(*it);
			}
			else {
				hestElem.dir[d] = &(*it);
			}
		}

		int numNodesAtFace;
		sm.Dimension() == 2 ? numNodesAtFace = numNodesAtFace = fec->GetOrder() + 1 : numNodesAtFace = getFaceNodeNumByGeomType(subFES);
		initNormalVectors(hestElem, numFaces * numNodesAtFace);
		initFscale(hestElem, numFaces * numNodesAtFace);

		auto elementVolume = sm.GetElementVolume(0);

		for (auto f{ 0 }; f < numFaces; f++) {
			
			Vector normal(sm.Dimension());
			ElementTransformation* faceTrans;
			sm.Dimension() == 2 ? faceTrans = sm.GetEdgeTransformation(f) : faceTrans = sm.GetFaceTransformation(f);
			faceTrans->SetIntPoint(&Geometries.GetCenter(faceTrans->GetGeometryType()));
			CalcOrtho(faceTrans->Jacobian(), normal);
			normal /= faceTrans->Weight();

			for (auto b{ 0 }; b < numNodesAtFace; b++) { //hesthaven requires normals to be stored once per node at face
				hestElem.normals[X][f * numNodesAtFace + b] = normal[0];
				hestElem.fscale[f * numNodesAtFace + b] = abs(normal[0] / elementVolume); //likewise for fscale, surface per volume ratio per node at face
				if (sm.Dimension() >= 2) {
					hestElem.normals[Y][f * numNodesAtFace + b] = normal[1];
					hestElem.fscale[f * numNodesAtFace + b] += abs(normal[1] / elementVolume);
				}
				if (sm.Dimension() == 3) {
					hestElem.normals[Z][f * numNodesAtFace + b] = normal[2];
					hestElem.fscale[f * numNodesAtFace + b] += abs(normal[2] / elementVolume);
				}
			}
		}

		hestElemStorage_[e] = hestElem;

	}

}

void HesthavenEvolution::Mult(const Vector& in, Vector& out) const
{
	double alpha;
	opts_.fluxType == FluxType::Upwind ? alpha = 1.0 : alpha = 0.0;

	// --MAP BETWEEN MFEM VECTOR AND EIGEN VECTOR-- //

	FieldsInputMaps fieldsIn(in, fes_);
	std::array<GridFunction, 3> eOut, hOut;

	for (int d = X; d <= Z; d++) {
		eOut[d].SetSpace(&fes_);
		hOut[d].SetSpace(&fes_);
		eOut[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hOut[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
		eOut[d] = 0.0;
		hOut[d] = 0.0;
	}

	// ---JUMPS--- //

	auto jumps{ HesthavenFields(connectivity_.size()) };

	for (auto v{ 0 }; v < connectivity_.size(); v++) {
		for (int d = X; d <= Z; d++) {
			jumps.e_[d][v] = fieldsIn.e_[d][connectivity_[v].first] - fieldsIn.e_[d][connectivity_[v].second];
			jumps.h_[d][v] = fieldsIn.h_[d][connectivity_[v].first] - fieldsIn.h_[d][connectivity_[v].second];
		}
	}

	// --BOUNDARIES-- //

	applyBoundaryConditionsToNodes(connectivity_, bdr_connectivity_, fieldsIn, jumps);

	// --INTERIOR BOUNDARIES-- //

	//if (int_bdr_connectivity_.size() != 0) {
	//	applyInteriorBoundaryConditionsToNodes(connectivity_, int_bdr_connectivity_ , fieldsIn, jumps);
	//}

	// --TOTAL FIELD SCATTERED FIELD-- //

	if (tfsf_connectivity_.size() != 0) {
		evaluateTFSF(jumps);
	}

	// --ELEMENT BY ELEMENT EVOLUTION-- //

	for (auto e{ 0 }; e < fes_.GetNE(); e++) {

		Array<int> dofs;
		auto el2dofs = fes_.GetElementDofs(e, dofs);
		auto elemFluxSize{ hestElemStorage_[e].fscale.size() };

		// Dof ordering will always be incremental due to L2 space (i.e: element 0 will have 0, 1, 2... element 1 will have 3, 4, 5...)

		const auto jumpsElem{ HesthavenElementJumps(jumps, e, elemFluxSize) };
		const auto fieldsElem{ FieldsElementMaps(in, fes_, e) };

		Eigen::VectorXd ndotdH(jumpsElem.h_[X].size()), ndotdE(jumpsElem.e_[X].size());
		
		for (int d = X; d <= Z; d++) {
			ndotdH += hestElemStorage_[e].normals[d].asDiagonal() * jumpsElem.h_[d];
			ndotdH += hestElemStorage_[e].normals[d].asDiagonal() * jumpsElem.e_[d];
		}

		HesthavenFields elemFlux(elemFluxSize);

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			elemFlux.h_[x] = -1.0 * hestElemStorage_[e].normals[y].asDiagonal() * jumpsElem.e_[z] +       hestElemStorage_[e].normals[z].asDiagonal() * jumpsElem.e_[y] + alpha * (jumpsElem.h_[x] - ndotdH.asDiagonal() * hestElemStorage_[e].normals[x]);
			elemFlux.e_[x] =        hestElemStorage_[e].normals[y].asDiagonal() * jumpsElem.h_[z] - 1.0 * hestElemStorage_[e].normals[z].asDiagonal() * jumpsElem.h_[y] + alpha * (jumpsElem.e_[x] - ndotdE.asDiagonal() * hestElemStorage_[e].normals[x]);

		}

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const auto& invmass = *hestElemStorage_[e].invmass;
			const auto& dir1 = *hestElemStorage_[e].dir[y];
			const auto& dir2 = *hestElemStorage_[e].dir[z];

			const DynamicMatrix result = invmass * dir1;
			
			const Eigen::VectorXd& hResult = -1.0 * invmass * dir1 * fieldsElem.e_[z] +       invmass * dir2 * fieldsElem.e_[y] + applyScalingFactors(hestElemStorage_[e], elemFlux.h_[x]);
			const Eigen::VectorXd& eResult =        invmass * dir1 * fieldsElem.h_[z] - 1.0 * invmass * dir2 * fieldsElem.h_[y] + applyScalingFactors(hestElemStorage_[e], elemFlux.e_[x]);
			const int startIndex = e * dofs.Size();
			mfem::real_t* mfemHFieldVals = new mfem::real_t[hResult.size()];
			mfem::real_t* mfemEFieldVals = new mfem::real_t[eResult.size()];
			for (auto v{ 0 }; v < hResult.size(); v++) {
				mfemHFieldVals[v] = hResult.data()[v];
				mfemEFieldVals[v] = eResult.data()[v];
			}
			hOut[x].SetSubVector(dofs, mfemHFieldVals);
			eOut[x].SetSubVector(dofs, mfemEFieldVals);

		}

	}

}

}

