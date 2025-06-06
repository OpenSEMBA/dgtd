#include "HesthavenEvolution.h"

#include <chrono>
#include <iostream>
#include <string>

namespace maxwell {

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
	const auto& mapBSF = connectivity_.boundary.TFSF.mapBSF;
	const auto& mapBTF = connectivity_.boundary.TFSF.mapBTF;
	const auto& vmapBSF = connectivity_.boundary.TFSF.vmapBSF;
	const auto& vmapBTF = connectivity_.boundary.TFSF.vmapBTF;
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
					out.e_[d][mapBSF[m][v]] -= fields[E][d];
					out.h_[d][mapBSF[m][v]] -= fields[H][d];
					out.e_[d][mapBTF[m][v]] += fields[E][d];
					out.h_[d][mapBTF[m][v]] += fields[H][d];
				}
			}
		}
	}

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

void HesthavenEvolution::storeDirectionalMatrices(ParFiniteElementSpace& subFES, const DynamicMatrix& refInvMass, HesthavenElement& hestElem)
{
	Model model(*subFES.GetMesh(), GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo{});
	Probes probes;
	ProblemDescription pd(model, probes, srcmngr_.sources, opts_);
	DGOperatorFactory dgops(pd, subFES);
	for (int d = X; d <= Z; d++) {
		auto denseMat = dgops.buildDerivativeSubOperator(d)->SpMat().ToDenseMatrix();
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

void HesthavenEvolution::applyBoundaryConditionsToNodes(const BoundaryMaps& bdrMaps, const FieldsInputMaps& in, HesthavenFields& out) const
{
	for (int m = 0; m < bdrMaps.PEC.vmapB.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.PEC.vmapB[m].size(); v++) {
				if (!isDoFinCurvedElement(bdrMaps.PEC.vmapB[m][v])) {
					out.e_[d][bdrMaps.PEC.mapB[m][v]] = -2.0 * in.e_[d][bdrMaps.PEC.vmapB[m][v]];
					out.h_[d][bdrMaps.PEC.mapB[m][v]] = 0.0;
				}
			}
		}
	}

	for (int m = 0; m < bdrMaps.PMC.vmapB.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.PMC.vmapB[m].size(); v++) {
				if (!isDoFinCurvedElement(bdrMaps.PMC.vmapB[m][v])) {
					out.e_[d][bdrMaps.PMC.mapB[m][v]] = 0.0;
					out.h_[d][bdrMaps.PMC.mapB[m][v]] = -2.0 * in.h_[d][bdrMaps.PMC.vmapB[m][v]];
				}
			}
		}
	}

	for (int m = 0; m < bdrMaps.SMA.mapB.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.SMA.vmapB[m].size(); v++) {
				if (!isDoFinCurvedElement(bdrMaps.SMA.vmapB[m][v])) {
					out.e_[d][bdrMaps.SMA.mapB[m][v]] = -1.0 * in.e_[d][bdrMaps.SMA.vmapB[m][v]] / opts_.alpha;
					out.h_[d][bdrMaps.SMA.mapB[m][v]] = -1.0 * in.h_[d][bdrMaps.SMA.vmapB[m][v]] / opts_.alpha;
				}
			}
		}
	}

	for (int m = 0; m < bdrMaps.intPEC.mapBElem1.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.intPEC.mapBElem1[m].size(); v++) { //Condition is applied twice, so we halve the coefficients for interior operators
				if (!isDoFinCurvedElement(bdrMaps.intPEC.vmapBElem1[m][v]) || !isDoFinCurvedElement(bdrMaps.intPEC.vmapBElem2[m][v])) {
					out.e_[d][bdrMaps.intPEC.mapBElem1[m][v]] = -1.0 * in.e_[d][bdrMaps.intPEC.vmapBElem1[m][v]];
					out.e_[d][bdrMaps.intPEC.mapBElem2[m][v]] = -1.0 * in.e_[d][bdrMaps.intPEC.vmapBElem2[m][v]];
					out.h_[d][bdrMaps.intPEC.mapBElem1[m][v]] = 0.0;
					out.h_[d][bdrMaps.intPEC.mapBElem2[m][v]] = 0.0;
				}
			}
		}
	}

	for (int m = 0; m < bdrMaps.intPMC.mapBElem1.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.intPMC.mapBElem1[m].size(); v++) {
				if (!isDoFinCurvedElement(bdrMaps.intPMC.vmapBElem1[m][v]) || !isDoFinCurvedElement(bdrMaps.intPMC.vmapBElem2[m][v])) {
					out.e_[d][bdrMaps.intPMC.mapBElem1[m][v]] = 0.0;
					out.e_[d][bdrMaps.intPMC.mapBElem2[m][v]] = 0.0;
					out.h_[d][bdrMaps.intPMC.mapBElem1[m][v]] = -1.0 * in.h_[d][bdrMaps.intPMC.vmapBElem1[m][v]];
					out.h_[d][bdrMaps.intPMC.mapBElem2[m][v]] = -1.0 * in.h_[d][bdrMaps.intPMC.vmapBElem2[m][v]];
				}
			}
		}
	}

	for (int m = 0; m < bdrMaps.intSMA.mapBElem1.size(); m++) {
		for (int d = X; d <= Z; d++) {
			for (int v = 0; v < bdrMaps.intSMA.mapBElem1[m].size(); v++) {
				if (!isDoFinCurvedElement(bdrMaps.intSMA.vmapBElem1[m][v]) || !isDoFinCurvedElement(bdrMaps.intSMA.vmapBElem2[m][v])) {
					out.e_[d][bdrMaps.intSMA.mapBElem1[m][v]] = -0.5 * in.e_[d][bdrMaps.intSMA.vmapBElem1[m][v]];
					out.h_[d][bdrMaps.intSMA.mapBElem1[m][v]] = -0.5 * in.h_[d][bdrMaps.intSMA.vmapBElem1[m][v]];
					out.e_[d][bdrMaps.intSMA.mapBElem2[m][v]] = -0.5 * in.e_[d][bdrMaps.intSMA.vmapBElem2[m][v]];
					out.h_[d][bdrMaps.intSMA.mapBElem2[m][v]] = -0.5 * in.h_[d][bdrMaps.intSMA.vmapBElem2[m][v]];
				}
			}
		}
	}

}

HesthavenEvolution::HesthavenEvolution(ParFiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, EvolutionOptions& opts) :
	TimeDependentOperator(numberOfFieldComponents* numberOfMaxDimensions* fes.GetNDofs()),
	fes_(fes),
	model_(model),
	srcmngr_(srcmngr),
	opts_(opts),
	connectivity_(model, fes)
{

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Hesthaven Evolution Operator is being initialized." << std::endl;
	std::cout << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
#endif

	Array<int> elementMarker;
	elementMarker.Append(hesthavenMeshingTag);

	const auto* cmesh = &model_.getConstMesh();
	auto mesh{ ParMesh(model_.getMesh()) };
	auto fec{ dynamic_cast<const L2_FECollection*>(fes_.FEColl()) };
	auto attMap{ mapOriginalAttributes(model_.getMesh()) };

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

#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << std::to_string(linearElements_.Size()) + " linear elements out of " + std::to_string(cmesh->GetNE()) + " total elements." << std::endl;
	std::cout << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
#endif

	for (int e = 0; e < linearElements_.Size(); e++)
	{
		HesthavenElement hestElem;
		hestElem.id = linearElements_[e];
		hestElem.type = cmesh->GetElementType(linearElements_[e]);
		hestElem.vol = mesh.GetElementVolume(linearElements_[e]);

		mesh.SetAttribute(linearElements_[e], hesthavenMeshingTag);
		auto sm{ ParSubMesh::CreateFromDomain(mesh, elementMarker) };
		restoreOriginalAttributesAfterSubMeshing(linearElements_[e], mesh, attMap);
		ParFiniteElementSpace subFES(&sm, fec);

		sm.bdr_attributes.SetSize(subFES.GetNF());
		for (auto f= 0; f < subFES.GetNF(); f++) {
			sm.bdr_attributes[f] = f + 1;
			sm.SetBdrAttribute(f, sm.bdr_attributes[f]);
		}

		storeDirectionalMatrices(subFES, refInvMass, hestElem);

		storeFaceInformation(subFES, hestElem);

		hestElemLinearStorage_[e] = hestElem;
	}
	
#ifdef SHOW_TIMER_INFORMATION
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << std::to_string(curvedElements_.size()) + " linear elements out of " + std::to_string(cmesh->GetNE()) + " total elements." << std::endl;
	std::cout << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
#endif

	if (curvedElements_.size()) {
		Probes probes;
		ProblemDescription pd(model_, probes, srcmngr_.sources, opts_);
		DGOperatorFactory dgops(pd, fes_);
		auto global = dgops.buildGlobalOperator();

		for (const auto& [e, dofs]: curvedElements_) {
			HesthavenCurvedElement hestCurElem;
			hestCurElem.id = e;
			hestCurElem.type = cmesh->GetElementType(e);
			hestCurElem.dofs = dofs;
			SparseMatrix spmat = SparseMatrix(numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs(), numberOfFieldComponents * numberOfMaxDimensions * fes_.GetNDofs());
			Array<int> cols; 
			Vector vals;
			for (auto d : dofs) {
				for (auto ft= 0; ft < numberOfFieldComponents * numberOfMaxDimensions; ft++) {
					global->GetRow(d + ft * fes_.GetNDofs(), cols, vals);
					spmat.SetRow(d + ft * fes_.GetNDofs(), cols, vals);
				}
			}
			spmat.Finalize();
			hestCurElem.matrix = std::move(spmat);
			hestElemCurvedStorage_.push_back(hestCurElem);
		}
	}

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

void HesthavenEvolution::Mult(const Vector& in, Vector& out) const
{
	double alpha = opts_.alpha;
	in.UseDevice(true);
	out.UseDevice(true);

	// --MAP BETWEEN MFEM VECTOR AND EIGEN VECTOR-- //

	const FieldsInputMaps fieldsIn(in, fes_);
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

	for (int v = 0; v < connectivity_.global.size(); v++) {
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

	for (int e = 0 ; e < linearElements_.Size(); e++) {
		auto elemFluxSize{ hestElemLinearStorage_[e].fscale.size() };

		// Dof ordering will always be incremental due to L2 space (i.e: element 0 will have 0, 1, 2... element 1 will have 3, 4, 5...)

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

			loadOutVectors(hResult, fes_, hestElemLinearStorage_[e].id, hOut[x]);
			loadOutVectors(eResult, fes_, hestElemLinearStorage_[e].id, eOut[x]);
		}

	}

	for (int e = 0; e < hestElemCurvedStorage_.size(); e++) {
		hestElemCurvedStorage_[e].matrix.AddMult(in, out);
    }

}
}