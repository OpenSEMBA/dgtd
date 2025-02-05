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

Mesh getMeshForReferenceElementBasedOnGeomType(const Element::Type elType, const int dimension)
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
		return Mesh::MakeCartesian3D(1, 1, 1, elType);
	default:
		throw std::runtime_error("Hesthaven Evolution Operator only supports dimensions 2 or 3.");
	}
}

DynamicMatrix assembleHesthavenReferenceElementInverseMassMatrix(const Element::Type elType, const int order, const int dimension)
{
	auto m{ getMeshForReferenceElementBasedOnGeomType(elType, dimension) };
	auto fec{ L2_FECollection(order, dimension, BasisType::GaussLobatto) };
	auto fes{ FiniteElementSpace(&m, &fec)};
	auto mass_mat{ assembleInverseMassMatrix(fes) };

	return getElementMassMatrixFromGlobal(0, mass_mat);
}

DynamicMatrix assembleHesthavenReferenceElementEmat(const Element::Type elType, const int order, const int dimension)
{
	auto m{ getMeshForReferenceElementBasedOnGeomType(elType, dimension) };
	m.SetAttribute(0, hesthavenMeshingTag);
	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);
	auto sm{ SubMesh::CreateFromDomain(m, elementMarker) };
	auto fec{ L2_FECollection(order, dimension, BasisType::GaussLobatto) };
	FiniteElementSpace subFES(&sm, &fec);

	auto boundary_markers = assembleBoundaryMarkers(subFES);

	for (auto f{ 0 }; f < subFES.GetNF(); f++) {
		sm.bdr_attributes[f] = f + 1;
		sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
	}

	return assembleEmat(subFES, boundary_markers);
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

void HesthavenEvolution::assembleTFSFConnectivity(const BilinearForm* matrix, FaceElementTransformations* faceTrans, double faceOri)
{
	Array<int> cols;
	Vector vals;
	double tol{ 1e-6 };
	for (auto r{ 0 }; r < matrix->NumRows() / 2; r++) {
		matrix->SpMat().GetRow(r, cols, vals);
		Nodes vmapBelem1, vmapBelem2;
		for (auto c{ 0 }; c < cols.Size(); c++) {
			if (std::abs(vals[c]) > tol && cols[c] != matrix->NumCols() - 1) {
				for (auto c2{ c + 1 }; c2 < cols.Size(); c2++) {
					if (std::abs(vals[c2]) > tol) {
						vmapBelem1.push_back((faceTrans->Elem1No * int(matrix->NumCols()) / 2) + cols[c]);
						vmapBelem2.push_back((faceTrans->Elem2No * int(matrix->NumCols()) / 2) + cols[c2] - (int(matrix->NumCols()) / 2));
					}
				}
			}
		}
		faceOri >= 0.0 ? bdr_connectivity_.TFSF.signs.push_back(std::make_pair(-1.0, 1.0)) : bdr_connectivity_.TFSF.signs.push_back(std::make_pair(-1.0, 1.0)); // Face Ori >= 0.0 means elem1 is SF, elem2 is TF.
		Nodes mapBElem1(vmapBelem1.size()), mapBElem2(vmapBelem2.size());
		for (auto v{ 0 }; v < mapBElem1.size(); v++) {
			mapBElem1[v] = std::distance(std::begin(connectivity_), std::find(connectivity_.begin(), connectivity_.end(), std::make_pair(vmapBelem1[v], vmapBelem2[v])));
			mapBElem2[v] = std::distance(std::begin(connectivity_), std::find(connectivity_.begin(), connectivity_.end(), std::make_pair(vmapBelem2[v], vmapBelem1[v])));
		}
		bdr_connectivity_.TFSF.mapBElem1.push_back(mapBElem1);
		bdr_connectivity_.TFSF.mapBElem2.push_back(mapBElem2);
		bdr_connectivity_.TFSF.vmapBElem1.push_back(vmapBelem1);
		bdr_connectivity_.TFSF.vmapBElem2.push_back(vmapBelem2);
	}
}

void HesthavenEvolution::evaluateTFSF(HesthavenFields& jumps) const
{
	auto func{ evalTimeVarFunction(GetTime(),srcmngr_) };
	for (auto m{ 0 }; m < bdr_connectivity_.TFSF.mapBElem1.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (auto v{ 0 }; v < bdr_connectivity_.TFSF.mapBElem1.size(); v++) {
				jumps.e_[d][bdr_connectivity_.TFSF.mapBElem1[m][v]] += func[E][d][bdr_connectivity_.TFSF.vmapBElem1[m][v]] * bdr_connectivity_.TFSF.signs[m].first;
				jumps.h_[d][bdr_connectivity_.TFSF.mapBElem1[m][v]] += func[H][d][bdr_connectivity_.TFSF.vmapBElem1[m][v]] * bdr_connectivity_.TFSF.signs[m].first;
				jumps.e_[d][bdr_connectivity_.TFSF.mapBElem2[m][v]] += func[E][d][bdr_connectivity_.TFSF.vmapBElem2[m][v]] * bdr_connectivity_.TFSF.signs[m].second;
				jumps.h_[d][bdr_connectivity_.TFSF.mapBElem2[m][v]] += func[H][d][bdr_connectivity_.TFSF.vmapBElem2[m][v]] * bdr_connectivity_.TFSF.signs[m].second;
			}
		}
	}
}

const Eigen::VectorXd HesthavenEvolution::applyScalingFactors(const ElementId e, const Eigen::VectorXd& flux) const
{
	const auto& fscale = hestElemStorage_[e].fscale.asDiagonal();
	return this->refLIFT_ * (fscale * flux) / 2.0;
}

std::vector<Array<int>> assembleBdrMarkersForBdrElements(Mesh& bdrMesh, const int numBdrElems)
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

Nodes findNodesPerBdrFace(const BilinearForm* bdrNodeMat)
{
	Nodes res;
	Array<int> cols;
	Vector vals;
	double tol{ 1e-6 };
	for (auto r{ 0 }; r < bdrNodeMat->NumRows(); r++) {
		bdrNodeMat->SpMat().GetRow(r, cols, vals);
		if (vals.Size() && std::abs(vals[cols.Find(r)]) > tol) {
			for (auto c{ 0 }; c < cols.Size(); c++) {
				if (std::abs(vals[c]) > 1e-5) {
					res.push_back(cols[c]);
				}
			}
			break;
		}
	}
	return res;
}

std::vector<Nodes> assembleNodeVectorPerBdrFace(std::vector<Array<int>>& bdrNodeMarkers, FiniteElementSpace& bdrFES, const std::map<bool, std::vector<BdrElementId>>& isInteriorMap)
{
	
	std::vector<Nodes> res(isInteriorMap.at(false).size());

	for (auto b{ 0 }; b < isInteriorMap.at(false).size(); b++) {
		auto bdrNodeFinderOperator{ assembleFaceMassBilinearForm(bdrFES, bdrNodeMarkers[isInteriorMap.at(false)[b]]) };
		res[b] = findNodesPerBdrFace(bdrNodeFinderOperator.get());
	}

	return res;
}

Nodes initMapB(const GlobalConnectivity& connectivity)
{
	Nodes res;
	for (auto i{ 0 }; i < connectivity.size(); i++) {
		if (connectivity[i].first == connectivity[i].second) {
			res.push_back(i);
		}
	}
	return res;
}


Nodes initVMapB(const GlobalConnectivity& connectivity, const Nodes& mapB)
{
	Nodes res(mapB.size());
	for (auto i{ 0 }; i < mapB.size(); i++) {
		res[i] = connectivity[mapB[i]].first;
	}
	return res;
}

const std::map<bool, std::vector<BdrElementId>> assembleInteriorOrTrueBdrMap(const FiniteElementSpace& fes)
{
	std::map<bool, std::vector<BdrElementId>> res;
	auto f2bdr{ fes.GetMesh()->GetFaceToBdrElMap() };
	for (auto b{ 0 }; b < fes.GetNBE(); b++) {
		fes.GetMesh()->FaceIsInterior(f2bdr.Find(b)) ==
			true ? res[true].push_back(b) : res[false].push_back(b);
	}
	return res;
}

void HesthavenEvolution::initBdrConnectivityMaps(const std::vector<Nodes>& bdr2nodes)
{
	auto mapB{ initMapB(connectivity_) };
	auto vmapB{ initVMapB(connectivity_, mapB) };

	auto modelBdrMarkers = model_.getBoundaryToMarker(); //Only for true boundary types (not interior)...
	for (auto b{ 0 }; b < bdr2nodes.size(); b++) { //For each one of our bdr elements...
		for (const auto& marker : modelBdrMarkers) { //Fetch each one of our bdrMarkers in the model...
			if (marker.second[fes_.GetBdrAttribute(b) - 1] == 1) { //And check if that BdrCond is active for the specified bdrAtt of that bdr element, which means the bdr element is of that type...
				for (auto n{ 0 }; n < vmapB.size(); n++) { //Then sweep vmapB...
					if (bdr2nodes[b][0] == connectivity_[mapB[n]].first) { //To find the nodes in vmapB that are equal to the first node in our bdr element (as they are sorted when built), because it can happen...
						Nodes mapBToStore(bdr2nodes[b].size()); //... that we have multiple instances of a Hesthaven type node appearing twice on vmapB (i.e. corner node)...
						Nodes vmapBToStore(bdr2nodes[b].size()); //... so we need to save the mapB 2 vmapB pairs so the nodes are properly linked when defining boundary conditions...
						for (auto m{ 0 }; m < bdr2nodes[b].size(); m++) {
							mapBToStore[m] = mapB[n - m]; //Then make a temporary vector that we'll fill the next nodesize worth of nodes for mapB and vmapB, starting at the compared and equal node...
							vmapBToStore[m] = connectivity_[mapB[n - m]].first;
						}
						if (bdr2nodes[b] == vmapBToStore) { //As bdr2nodes[b] is already sorted when built, potVec also has to be sorted to match...
							switch (marker.first) {
							case BdrCond::PEC:
								bdr_connectivity_.PEC.mapB.push_back(mapBToStore);
								bdr_connectivity_.PEC.vmapB.push_back(vmapBToStore);
								break;
							case BdrCond::PMC:
								bdr_connectivity_.PMC.mapB.push_back(mapBToStore);
								bdr_connectivity_.PMC.vmapB.push_back(vmapBToStore);
								break;
							case BdrCond::SMA:
								bdr_connectivity_.SMA.mapB.push_back(mapBToStore);
								bdr_connectivity_.SMA.vmapB.push_back(vmapBToStore);
								break;
							}
							break; // And as we only support one bdrcond per face/edge/bdrwhatever, we're done for this bdr element...
						}     // As this method is completely non-dependent on dimensions or geometries, could be assumed is as generic as it could get.
					}
				}
			}
		}
	}
}

void assembleDerivativeMatrices(FiniteElementSpace& subFES, MatrixStorageLT& matStLt, HesthavenElement& hestElem)
{
	for (auto d{ X }; d <= Z; d++) {
		auto derivativeMatrix{ toEigen(*buildDerivativeOperator(d, subFES)->SpMat().ToDenseMatrix()) };
		StorageIterator it = matStLt.find(derivativeMatrix);
		if (it == matStLt.end()) {
			matStLt.insert(derivativeMatrix);
			StorageIterator it = matStLt.find(derivativeMatrix);
			hestElem.der[d] = &(*it);
		}
		else {
			hestElem.der[d] = &(*it);
		}
	}
}

void assembleFaceInformation(FiniteElementSpace& subFES, HesthavenElement& hestElem)
{

	int numFaces, numNodesAtFace;
	const auto& dim = subFES.GetMesh()->Dimension();
	dim == 2 ? numFaces = subFES.GetMesh()->GetNEdges() : numFaces = subFES.GetMesh()->GetNFaces();
	dim == 2 ? numNodesAtFace = numNodesAtFace = subFES.FEColl()->GetOrder() + 1 : numNodesAtFace = getFaceNodeNumByGeomType(subFES);
	initNormalVectors(hestElem, numFaces * numNodesAtFace);
	initFscale(hestElem, numFaces * numNodesAtFace);

	for (auto f{ 0 }; f < numFaces; f++) {

		Vector normal(dim);
		ElementTransformation* faceTrans;
		dim == 2 ? faceTrans = subFES.GetMesh()->GetEdgeTransformation(f) : faceTrans = subFES.GetMesh()->GetFaceTransformation(f);
		faceTrans->SetIntPoint(&Geometries.GetCenter(faceTrans->GetGeometryType()));
		CalcOrtho(faceTrans->Jacobian(), normal);
		const auto sJ{ faceTrans->Weight() };

		for (auto b{ 0 }; b < numNodesAtFace; b++) { //hesthaven requires normals to be stored once per node at face
			hestElem.normals[X][f * numNodesAtFace + b] = normal[0] / sJ;
			hestElem.fscale[f * numNodesAtFace + b] = sJ / hestElem.vol; //likewise for fscale, surface per volume ratio per node at face
			if (dim >= 2) {
				hestElem.normals[Y][f * numNodesAtFace + b] = normal[1] / sJ;
			}
			if (dim == 3) {
				hestElem.normals[Z][f * numNodesAtFace + b] = normal[2] / sJ;
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

	const auto* cmesh = &model_.getConstMesh();
	auto mesh{ Mesh(model_.getMesh()) };
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
	auto attMap{ mapOriginalAttributes(model_.getMesh()) };

	const auto isInteriorMap{ assembleInteriorOrTrueBdrMap(fes_) };
	
	connectivity_ =  assembleGlobalConnectivityMap(mesh, fec);

	{
		auto bdrNodeMesh{ Mesh(model_.getMesh()) };
		std::vector<Array<int>> bdrNodeMarkers{ assembleBdrMarkersForBdrElements(bdrNodeMesh, model_.getConstMesh().GetNBE()) };
		auto bdrNodeFES = FiniteElementSpace(&bdrNodeMesh, fec);
		const auto bdr2nodes{ assembleNodeVectorPerBdrFace(bdrNodeMarkers, bdrNodeFES, isInteriorMap) };
		initBdrConnectivityMaps(bdr2nodes);
	}

	if (model_.getInteriorBoundaryToMarker().size() != 0) {
		auto intBdrNodeMesh{ Mesh(model_.getMesh()) };
		for (auto b{ 0 }; b < model_.getConstMesh().GetNBE(); b++) {
			for (const auto& marker : model_.getInteriorBoundaryToMarker()) {
				if (marker.second[model_.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
					auto faceTrans{ model_.getMesh().GetInternalBdrFaceTransformations(b) };
					intBdrNodeMesh.SetAttribute(faceTrans->Elem1No, hesthavenMeshingTag);
					intBdrNodeMesh.SetAttribute(faceTrans->Elem2No, hesthavenMeshingTag);
					auto twoElementSubMesh = SubMesh::CreateFromDomain(intBdrNodeMesh, elementMarker);
					restoreOriginalAttributesAfterSubMeshing(faceTrans->Elem1No, intBdrNodeMesh, attMap);
					restoreOriginalAttributesAfterSubMeshing(faceTrans->Elem2No, intBdrNodeMesh, attMap);
					FiniteElementSpace subFES(&twoElementSubMesh, fec);
					auto nodePairs{ buildConnectivityForInteriorBdrFace(*faceTrans, fes_, subFES) };
					Nodes mapBElem1(nodePairs.first.size()), mapBElem2(nodePairs.first.size());
					for (auto v{ 0 }; v < mapBElem1.size(); v++) {
						mapBElem1[v] = std::distance(std::begin(connectivity_), std::find(connectivity_.begin(), connectivity_.end(), std::make_pair(nodePairs.first[v], nodePairs.second[v])));
						mapBElem2[v] = std::distance(std::begin(connectivity_), std::find(connectivity_.begin(), connectivity_.end(), std::make_pair(nodePairs.second[v], nodePairs.first[v])));
					}
					if (marker.first != BdrCond::TotalFieldIn) {
						switch (marker.first) {
							case BdrCond::PEC:
								bdr_connectivity_.intPEC.mapBElem1.push_back(mapBElem1);
								bdr_connectivity_.intPEC.mapBElem2.push_back(mapBElem2);
								bdr_connectivity_.intPEC.vmapBElem1.push_back(nodePairs.first);
								bdr_connectivity_.intPEC.vmapBElem2.push_back(nodePairs.second);
								break;
							case BdrCond::PMC:
								bdr_connectivity_.intPMC.mapBElem1.push_back(mapBElem1);
								bdr_connectivity_.intPMC.mapBElem2.push_back(mapBElem2);
								bdr_connectivity_.intPMC.vmapBElem1.push_back(nodePairs.first);
								bdr_connectivity_.intPMC.vmapBElem2.push_back(nodePairs.second);
								break;
							case BdrCond::SMA:
								bdr_connectivity_.intSMA.mapBElem1.push_back(mapBElem1);
								bdr_connectivity_.intSMA.mapBElem2.push_back(mapBElem2);
								bdr_connectivity_.intSMA.vmapBElem1.push_back(nodePairs.first);
								bdr_connectivity_.intSMA.vmapBElem2.push_back(nodePairs.second);
								break;
							default:
								throw std::runtime_error("Incorrect boundary condition for interior assignment.");
						}
					}
				}
			}
		}
	}

	if (model_.getTotalFieldScatteredFieldToMarker().size() != 0) {
		auto mesh_tfsf{ Mesh(model_.getMesh()) };
		for (auto b{ 0 }; b < mesh_tfsf.GetNBE(); b++) {
			for (const auto& marker : model_.getTotalFieldScatteredFieldToMarker()) {
				if (marker.second[model_.getConstMesh().GetBdrAttribute(b) - 1] == 1) {
					auto faceOri{ buildFaceOrientation(mesh, b) };
					auto faceTrans{ mesh_tfsf.GetInternalBdrFaceTransformations(b) };
					auto sm{ assembleInteriorFaceSubMesh(mesh, *faceTrans, attMap) };
					FiniteElementSpace smFES(&sm, fec);
					auto matrix{ assembleInteriorFluxMatrix(smFES) };
					assembleTFSFConnectivity(matrix.get(), faceTrans, faceOri);
				}
			}
		}
	}

	mesh = Mesh(model_.getMesh()); 

	hestElemStorage_.resize(cmesh->GetNE());

	bool allElementsSameGeom = true;
	{
		const auto firstElemGeom = cmesh->GetElementGeometry(0);
		for (auto e{ 0 }; e < cmesh->GetNE(); e++)
		{
			if (firstElemGeom != cmesh->GetElementGeometry(e))
			{
				allElementsSameGeom = false;
			}
		}
	}
	
	if (allElementsSameGeom) 
	{
		refInvMass_ = assembleHesthavenReferenceElementInverseMassMatrix(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
		refLIFT_ = refInvMass_ * assembleHesthavenReferenceElementEmat(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
	}

	for (auto e{ 0 }; e < cmesh->GetNE(); e++)
	{
		HesthavenElement hestElem;
		hestElem.id = e;
		hestElem.geom = cmesh->GetElementBaseGeometry(e);
		hestElem.vol = mesh.GetElementVolume(e);

		mesh.SetAttribute(e, hesthavenMeshingTag);
		auto sm{ SubMesh::CreateFromDomain(mesh, elementMarker) };
		restoreOriginalAttributesAfterSubMeshing(e, mesh, attMap);
		FiniteElementSpace subFES(&sm, fec);

		sm.bdr_attributes.SetSize(subFES.GetNF());
		for (auto f{ 0 }; f < subFES.GetNF(); f++) {
			sm.bdr_attributes[f] = f + 1;
			sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
		}

		assembleDerivativeMatrices(subFES, matrixStorage_, hestElem);

		assembleFaceInformation(subFES, hestElem);

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
			jumps.e_[d][v] = fieldsIn.e_[d][connectivity_[v].second] - fieldsIn.e_[d][connectivity_[v].first];
			jumps.h_[d][v] = fieldsIn.h_[d][connectivity_[v].second] - fieldsIn.h_[d][connectivity_[v].first];
		}
	}

	// --BOUNDARIES-- //

	applyBoundaryConditionsToNodes(bdr_connectivity_, fieldsIn, jumps);

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
		ndotdH.setZero(); ndotdE.setZero();

		for (int d = X; d <= Z; d++) {
			ndotdH += hestElemStorage_[e].normals[d].asDiagonal() * jumpsElem.h_[d];
			ndotdE += hestElemStorage_[e].normals[d].asDiagonal() * jumpsElem.e_[d];
		}

		HesthavenFields elemFlux(elemFluxSize);

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const auto& norx = hestElemStorage_[e].normals[x];
			const auto& nory = hestElemStorage_[e].normals[y];
			const auto& norz = hestElemStorage_[e].normals[z];

			elemFlux.h_[x] = -1.0 * nory.asDiagonal() * jumpsElem.e_[z] +       norz.asDiagonal() * jumpsElem.e_[y] + alpha * (jumpsElem.h_[x] - ndotdH.asDiagonal() * norx);
			elemFlux.e_[x] =        nory.asDiagonal() * jumpsElem.h_[z] - 1.0 * norz.asDiagonal() * jumpsElem.h_[y] + alpha * (jumpsElem.e_[x] - ndotdE.asDiagonal() * norx);

		}

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const DynamicMatrix& invmass = this->refInvMass_ * (2.0 / hestElemStorage_[e].vol);
			const auto& der1 = *hestElemStorage_[e].der[y];
			const auto& der2 = *hestElemStorage_[e].der[z];

			const Eigen::VectorXd& hResult = -1.0 * invmass * der1 * fieldsElem.e_[z] +       invmass * der2 * fieldsElem.e_[y] + applyScalingFactors(e, elemFlux.h_[x]);
			const Eigen::VectorXd& eResult =        invmass * der1 * fieldsElem.h_[z] - 1.0 * invmass * der2 * fieldsElem.h_[y] + applyScalingFactors(e, elemFlux.e_[x]);

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

