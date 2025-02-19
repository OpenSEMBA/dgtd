#include "Evolution.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

using MatricesSet = std::set<DynamicMatrix, MatrixCompareLessThan>;


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
	auto fes{ FiniteElementSpace(&m, &fec)};
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
	for (auto f{ 0 }; f < subFES.GetNF(); f++) {
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
	const auto& mapBSF = connectivity_.boundary.TFSF.mapBSF;
	const auto& mapBTF = connectivity_.boundary.TFSF.mapBTF;
	const auto& vmapBSF = connectivity_.boundary.TFSF.vmapBSF;
	const auto& vmapBTF = connectivity_.boundary.TFSF.vmapBTF;
	for (const auto& source : srcmngr_.sources) {
		auto pw = dynamic_cast<Planewave*>(source.get());
		if (pw == nullptr) {
			continue;
		}
		for (auto v = 0; v < positions_.size(); v++) {
			if (std::abs(positions_[v][0] - 1.0) < 1e-2 && std::abs(positions_[v][1] - 0.0) < 1e-2) {
				int a = 0;
			}
		}
		for (auto m{ 0 }; m < mapBSF.size(); m++) {
			for (auto v{ 0 }; v < mapBSF[m].size(); v++){
				for (auto d : { X, Y, Z }) {
					fields[E][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), E, d);
					fields[H][d] = source->eval(positions_[vmapBSF[m][v]], GetTime(), H, d);
					out.e_[d][mapBSF[m][v]] -= fields[E][d];
					out.h_[d][mapBSF[m][v]] -= fields[H][d];
					out.e_[d][mapBTF[m][v]] += fields[E][d];
					out.h_[d][mapBTF[m][v]] += fields[H][d];
				}
			}
		}
	}

}

const Eigen::VectorXd HesthavenEvolution::applyLIFT(const ElementId e, Eigen::VectorXd& flux) const
{
	for (auto i{0}; i < flux.size(); i++) {
		flux[i] *= hestElemStorage_[e].fscale[i] / 2.0 ;
	}
	return this->refLIFT_ * flux;
}

void storeDirectionalMatrices(FiniteElementSpace& subFES, MatricesSet& matStLt, HesthavenElement& hestElem)
{
	for (auto d{ X }; d <= Z; d++) {
		auto denseMat = std::make_unique<mfem::DenseMatrix>(
			buildDerivativeOperator(d, subFES)->SpMat().ToDenseMatrix());
		auto derivativeMatrix{ toEigen(*denseMat) };
		StorageIterator it = matStLt.find(derivativeMatrix);
		if (it == matStLt.end()) {
			matStLt.insert(derivativeMatrix);
			StorageIterator it = matStLt.find(derivativeMatrix);
			hestElem.der[d] = &(*it);
		} else {
			hestElem.der[d] = &(*it);
		}
	}
}

double getReferenceVolume(const Element::Type geom)
{
	switch (geom) {
	case Element::Type::TRIANGLE:
		return 2.0; //Hesthaven definition (-1,-1), (-1,1), (1,-1)
	case Element::Type::QUADRILATERAL:
		return 4.0; //Assuming x,y (-1, 1)
	case Element::Type::TETRAHEDRON:
		return 8.0/6.0; //Hesthaven definition (-1,-1,-1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)
	case Element::Type::HEXAHEDRON:
		return 8.0; //Assuming x,y,z (-1, 1)
	default:
		throw std::runtime_error("Unsupported geometry for reference volume.");
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
	auto J{ subFES.GetMesh()->GetElementTransformation(0)->Weight()};

	for (auto f{ 0 }; f < numFaces; f++) {

		Vector normal(dim);
		ElementTransformation* faceTrans;
		dim == 2 ? faceTrans = subFES.GetMesh()->GetEdgeTransformation(f) : faceTrans = subFES.GetMesh()->GetFaceTransformation(f);
		CalcOrtho(faceTrans->Jacobian(), normal);
		auto sJ{ faceTrans->Weight() };

		for (auto b{ 0 }; b < numNodesAtFace; b++) { //hesthaven requires normals to be stored once per node at face
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

HesthavenEvolution::HesthavenEvolution(FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts),
	connectivity_(model, fes)
{
	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);

	const auto* cmesh = &model_.getConstMesh();
	auto mesh{ Mesh(model_.getMesh()) };
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
	auto attMap{ mapOriginalAttributes(model_.getMesh()) };	

	{
		if (model_.getTotalFieldScatteredFieldToMarker().size() != 0) {
			positions_ = buildDoFPositions(fes_);
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
		refInvMass_ = assembleHesthavenRefElemInvMassMatrix(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
		refLIFT_ = refInvMass_ * assembleHesthavenRefElemEmat(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
	}

	for (auto e{ 0 }; e < cmesh->GetNE(); e++)
	{
		HesthavenElement hestElem;
		hestElem.id = e;
		hestElem.type = cmesh->GetElementType(e);
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

		storeDirectionalMatrices(subFES, matrixStorage_, hestElem);

		assembleFaceInformation(subFES, hestElem);

		hestElemStorage_[e] = hestElem;

	}

}

void loadOutVectors(const Eigen::VectorXd& data, const FiniteElementSpace& fes, const ElementId& e, GridFunction& out) 
{
	Array<int> dofs;
	auto el2dofs = fes.GetElementDofs(e, dofs);
	std::unique_ptr<mfem::real_t[]> mfemFieldVars = std::make_unique<mfem::real_t[]>(data.size());
	for (auto v{ 0 }; v < data.size(); v++) {
		mfemFieldVars.get()[v] = data.data()[v];
	}
	out.SetSubVector(dofs, mfemFieldVars.get());
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

	auto jumps{ HesthavenFields(connectivity_.global.size()) };

	for (auto v{ 0 }; v < connectivity_.global.size(); v++) {
		for (int d = X; d <= Z; d++) {
			jumps.e_[d][v] = fieldsIn.e_[d][connectivity_.global[v].second] - fieldsIn.e_[d][connectivity_.global[v].first];
			jumps.h_[d][v] = fieldsIn.h_[d][connectivity_.global[v].second] - fieldsIn.h_[d][connectivity_.global[v].first];
		}
	}

	// --BOUNDARIES-- //

	applyBoundaryConditionsToNodes(connectivity_.boundary, fieldsIn, jumps);

	// --TOTAL FIELD SCATTERED FIELD-- //

	evaluateTFSF(jumps);

	// --ELEMENT BY ELEMENT EVOLUTION-- //

	for (auto e{ 0 }; e < fes_.GetNE(); e++) {

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

		auto refVol{ getReferenceVolume(hestElemStorage_[e].type) };
		const DynamicMatrix& invmass = this->refInvMass_ * refVol / hestElemStorage_[e].vol;

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const auto& der1 = *hestElemStorage_[e].der[y];
			const auto& der2 = *hestElemStorage_[e].der[z];

			const Eigen::VectorXd& hResult = -1.0 * invmass * der1 * fieldsElem.e_[z] +       invmass * der2 * fieldsElem.e_[y] + applyLIFT(e, elemFlux.h_[x]);
			const Eigen::VectorXd& eResult =        invmass * der1 * fieldsElem.h_[z] - 1.0 * invmass * der2 * fieldsElem.h_[y] + applyLIFT(e, elemFlux.e_[x]);

			loadOutVectors(hResult, fes_, e, hOut[x]);
			loadOutVectors(eResult, fes_, e, eOut[x]);
		}

	}

}

}

