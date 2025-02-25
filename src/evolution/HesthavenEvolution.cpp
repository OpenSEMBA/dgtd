#include "HesthavenEvolution.h"

namespace maxwell {

using MatricesSet = std::set<DynamicMatrix, MatrixCompareLessThan>;

static DynamicMatrix assembleInverseMassMatrix(FiniteElementSpace& fes)
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
			for (auto v{ 0 }; v < mapBSF[m].size(); v++) {
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
	for (auto i{ 0 }; i < flux.size(); i++) {
		flux[i] *= hestElemStorage_[e].fscale[i] / 2.0;
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
	Model model(*subFES.GetMesh(), GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo{});
	ProblemDescription pd(model, Probes{}, srcmngr_.sources, opts_);
	DGOperatorFactory dgops(pd, subFES);
	for (auto d{ X }; d <= Z; d++) {
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
	const auto& dim = subFES.GetMesh()->Dimension();
	dim == 2 ? numFaces = subFES.GetMesh()->GetNEdges() : numFaces = subFES.GetMesh()->GetNFaces();
	dim == 2 ? numNodesAtFace = numNodesAtFace = subFES.FEColl()->GetOrder() + 1 : numNodesAtFace = getFaceNodeNumByGeomType(subFES);
	initNormalVectors(hestElem, numFaces * numNodesAtFace);
	initFscale(hestElem, numFaces * numNodesAtFace);
	auto J{ subFES.GetMesh()->GetElementTransformation(0)->Weight() };

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

	DynamicMatrix refInvMass;
	if (allElementsSameGeom)
	{
		refInvMass = assembleHesthavenRefElemInvMassMatrix(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
		refLIFT_ = refInvMass * assembleHesthavenRefElemEmat(cmesh->GetElementType(0), fec->GetOrder(), cmesh->Dimension());
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

		storeDirectionalMatrices(subFES, refInvMass, hestElem);

		storeFaceInformation(subFES, hestElem);

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

			const Eigen::VectorXd& norx = hestElemStorage_[e].normals[x];
			const Eigen::VectorXd& nory = hestElemStorage_[e].normals[y];
			const Eigen::VectorXd& norz = hestElemStorage_[e].normals[z];

			elemFlux.h_[x] = -1.0 * nory.asDiagonal() * jumpsElem.e_[z] +       norz.asDiagonal() * jumpsElem.e_[y] + alpha * (jumpsElem.h_[x] - ndotdH.asDiagonal() * norx);
			elemFlux.e_[x] =        nory.asDiagonal() * jumpsElem.h_[z] - 1.0 * norz.asDiagonal() * jumpsElem.h_[y] + alpha * (jumpsElem.e_[x] - ndotdE.asDiagonal() * norx);

		}

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			const DynamicMatrix& dir1 = *hestElemStorage_[e].dir[y];
			const DynamicMatrix& dir2 = *hestElemStorage_[e].dir[z];

			const Eigen::VectorXd& hResult = -1.0 * dir1 * fieldsElem.e_[z] +        dir2 * fieldsElem.e_[y] + applyLIFT(e, elemFlux.h_[x]);
			const Eigen::VectorXd& eResult =        dir1 * fieldsElem.h_[z] - 1.0 *  dir2 * fieldsElem.h_[y] + applyLIFT(e, elemFlux.e_[x]);

			loadOutVectors(hResult, fes_, e, hOut[x]);
			loadOutVectors(eResult, fes_, e, eOut[x]);
		}

	}

}
}