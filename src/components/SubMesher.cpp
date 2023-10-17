#include "SubMesher.h"

namespace maxwell {

using namespace mfem;
static const int NotFound{ -1 };

TotalFieldScatteredFieldSubMesher::TotalFieldScatteredFieldSubMesher(const Mesh& m)
{
	Mesh parent_for_global(m);
	Mesh parent_for_individual(m);

	setGlobalTFSFAttributesForSubMeshing(parent_for_global);

	switch (m.Dimension()) {
	case 1:
		setIndividualTFSFAttributesForSubMeshing1D(parent_for_individual);
		break;
	default:
		setIndividualTFSFAttributesForSubMeshing(parent_for_individual);
		break;
	}

	Array<int> global_att(1); global_att[0] = SubMeshingMarkers::Global_SubMesh;
	auto global_sm{ SubMesh::CreateFromDomain(parent_for_global, global_att) };
	restoreElementAttributes(global_sm);
	global_sm.FinalizeMesh();
	global_submesh_ = std::make_unique<SubMesh>(global_sm);

	if (!elem_to_face_tf_.empty()) {
		tf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, true));
	}

	if (!elem_to_face_sf_.empty()) {
		sf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, false));
	}

};

SubMesh TotalFieldScatteredFieldSubMesher::createSubMeshFromParent(const Mesh& parent, bool isTF)
{
	Array<int> marker(1); 
	if (isTF) {
		marker[0] = SubMeshingMarkers::TotalField;
	}
	else {
		marker[0] = SubMeshingMarkers::ScatteredField;
	}
	auto res{ SubMesh::CreateFromDomain(parent, marker) };
	switch (parent.Dimension()) {
	case 1:
		setBoundaryAttributesInChild1D(parent, res);
		break;
	case 2:
		setBoundaryAttributesInChild2D(parent, res);
		break;
	case 3:
		setBoundaryAttributesInChild3D(parent, res);
		break;
	}
	restoreElementAttributes(res);
	res.FinalizeMesh();
	return res;
}

void TotalFieldScatteredFieldSubMesher::restoreElementAttributes(Mesh& m) //Temporary method that has to be reworked when materials are implemented.
{
	for (int e = 0; e < m.GetNE(); e++)
	{
		m.SetAttribute(e, 1);
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild1D(const Mesh& parent, SubMesh& child)
{
	for (int e = 0; e < child.GetNE(); e++) {
		Array<int> verts(2);
		child.GetElementVertices(e, verts);
		child.AddBdrPoint(verts[0]);
		child.AddBdrPoint(verts[1]);
	}
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int i = 0; i < parent.GetNBE(); i++) {
		if (parent.GetBdrAttribute(i) == static_cast<int>(BdrCond::TotalFieldIn)) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], static_cast<int>(BdrCond::TotalFieldIn));
		}
		else if (parent.GetBdrAttribute(i) == static_cast<int>(BdrCond::TotalFieldOut)) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], static_cast<int>(BdrCond::TotalFieldOut));
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild2D(const Mesh& parent, SubMesh& child)
{
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int i = 0; i < parent.GetNBE(); i++) {
		if (parent.GetBdrAttribute(i) == static_cast<int>(BdrCond::TotalFieldIn)) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(i))], static_cast<int>(BdrCond::TotalFieldIn));
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setBoundaryAttributesInChild3D(const Mesh& parent, SubMesh& child)
{
	setBoundaryAttributesInChild2D(parent, child);
}

void TotalFieldScatteredFieldSubMesher::setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf)
{
	if (trans->Elem2No != NotFound) {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalField);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::ScatteredField);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredField);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::TotalField);
		}
	}
	else {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalField);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredField);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf)
{
	if (facesId.second != NotFound) {
		if (el1_is_tf) {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
		}
		else {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem2No, facesId.second));
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
	}
	else {
		if (el1_is_tf) {
			elem_to_face_tf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
		else {
			elem_to_face_sf_.push_back(std::make_pair(trans->Elem1No, facesId.first));
		}
	}
}

void TotalFieldScatteredFieldSubMesher::prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf)
{
	setAttributeForTagging(m, trans, el1_is_tf);
	storeElementToFaceInformation(trans, facesId, el1_is_tf);
}

FaceElementTransformations* TotalFieldScatteredFieldSubMesher::getFaceElementTransformation(Mesh&m, int be) 
{
	switch (m.FaceIsInterior(m.GetBdrFace(be))) {
	case true:
		return m.GetInternalBdrFaceTransformations(be);
	default:
		return m.GetBdrFaceTransformations(be);
	}
}

void TotalFieldScatteredFieldSubMesher::setGlobalTFSFAttributesForSubMeshing(Mesh& m)
{

	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn) || 
			m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
			auto be_trans{ getFaceElementTransformation(m, be)};
			if (be_trans->Elem2No != NotFound) {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				m.GetElement(be_trans->Elem2No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
				elems_for_global_submesh_.push_back(be_trans->Elem2No);
			}
			else {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::Global_SubMesh);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
			}
		}
	}
}

Vector getBarycenterOfElement(Mesh& m, int e)
{
	Element* elem{ m.GetElement(e) };
	Array<int> elem_vert(elem->GetNVertices());
	elem->GetVertices(elem_vert);
	Vector res;
	switch (m.Dimension()) {
	case 1:
		res.SetSize(1);
		break;
	case 2:
		res.SetSize(2);
		break;
	case 3:
		res.SetSize(3);
		break;
	default:
		throw std::exception("Wrong dimension for TFSF barycenter calculation.");
	}
	res = 0.0;
	for (int v = 0; v < elem_vert.Size(); v++) {
		Vector vertexPos(m.GetVertex(elem_vert[v]), m.Dimension());
		res += vertexPos;
	}
	res /= elem_vert.Size();
	return res;
}

Vector subtract(const double* bdr_v, const Vector& b_v)
{
	Vector res(3);
	switch (b_v.Size()) {
	case 1:
		res.SetSize(1);
		break;
	case 2:
		res.SetSize(2);
		break;
	case 3:
		res.SetSize(3);
		break;
	default:
		throw std::exception("Wrong dimension for TFSF barycenter calculation.");
	}
	for (int i = 0; i < res.Size(); i++) {
		res[i] = b_v[i] - bdr_v[i];
	}
	return res;
}

Vector getNormal(FaceElementTransformations& fet)
{
	Vector res;
	switch (fet.Elem1->GetDimension()){
	case 1:
		res.SetSize(1);
		res(0) = 1.0;
		break;
	case 2: 
		res.SetSize(2);
		fet.SetIntPoint(&Geometries.GetCenter(fet.GetGeometryType()));
		CalcOrtho(fet.Jacobian(), res);
		break;
	case 3:
		res.SetSize(3);
		fet.SetIntPoint(&Geometries.GetCenter(fet.GetGeometryType()));
		CalcOrtho(fet.Jacobian(), res);
		break;
	default:
		throw std::exception("Wrong Dimension for Element in Normal for TFSF orientations.");
	}
	return res;
}

std::pair<double, double> calculateBaryNormalProduct(Mesh& m, FaceElementTransformations& fet, int be)
{
	auto bdr_vertices{ m.GetBdrElement(be)->GetVertices() };
	auto bdr_vert{ m.GetVertex(bdr_vertices[0]) };
	auto normal{ getNormal(fet) };
	
	auto bary1{ getBarycenterOfElement(m, fet.Elem1No) };
	auto v1{ subtract(bdr_vert,bary1) };
	auto d1{ mfem::InnerProduct(v1, normal) };

	auto d2 = 0.0;
	if (fet.Elem2No != NotFound) {
		auto bary2{ getBarycenterOfElement(m, fet.Elem2No) };
		auto v2{ subtract(bdr_vert,bary2) };
		d2 = mfem::InnerProduct(v2, normal);
	}

	return std::make_pair(d1, d2);
}

void TotalFieldScatteredFieldSubMesher::assignIndividualTFSFAttsOnePoint1D(Mesh& m)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {
			auto be_trans{ getFaceElementTransformation(m, be) };
			std::pair<FaceId, IsTF> set_e1;
			std::pair<FaceId, IsTF> set_e2;
			if (be_trans->Elem1No == 0) {
				set_e1 = std::make_pair(0, true);
				set_e2 = std::make_pair(NotFound, false);
			}
			else if (be_trans->Elem1No == m.GetNE() - 1) {
				set_e1 = std::make_pair(1, true);
				set_e2 = std::make_pair(NotFound, false);
			}
			else {
				set_e1 = std::make_pair(1, false);
				set_e2 = std::make_pair(0, true);
			}
			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_e1.first, set_e2.first);
			std::pair<IsTF, IsTF> dirInfo = std::make_pair(set_e1.second, set_e2.second);
			prepareSubMeshInfo(m, be_trans, facesInfo, set_e1.second);
		}
	}
}

SetPairs TotalFieldScatteredFieldSubMesher::twoPointAssignator(Mesh& m, int be, bool flag)
{
	auto be_trans{ getFaceElementTransformation(m, be) };
	std::pair<FaceId, IsTF> set_e1;
	std::pair<FaceId, IsTF> set_e2;
	switch (flag) {
	case false: //If first boundary we find
		if (be_trans->Elem1No == 0 && be_trans->Elem2No == NotFound) {
			set_e1 = std::make_pair(0, true);
			set_e2 = std::make_pair(NotFound, false);
		}
		else {
			set_e1 = std::make_pair(1, false);
			set_e2 = std::make_pair(0, true);
		}
		break;
	case true: //If second boundary we find
		if (be_trans->Elem1No == m.GetNE() - 1 && be_trans->Elem2No == NotFound) {
			set_e1 = std::make_pair(1, true);
			set_e2 = std::make_pair(NotFound, false);
		}
		else {
			set_e1 = std::make_pair(1, true);
			set_e2 = std::make_pair(0, false);
		}
		break;
	}
	return std::make_pair(set_e1, set_e2);
}

void TotalFieldScatteredFieldSubMesher::assignIndividualTFSFAttsTwoPoints1D(Mesh& m)
{
	auto flag{ false };
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {
			SetPairs sets;
			switch (flag) {
			case false:
				sets = twoPointAssignator(m, be, flag);
				flag = true;
				break;
			case true:
				sets = twoPointAssignator(m, be, flag);
				break;
			}

			std::pair<FaceId, FaceId> facesInfo = std::make_pair(sets.first.first, sets.second.first);
			std::pair<IsTF, IsTF> dirInfo = std::make_pair(sets.first.second, sets.second.second);
			prepareSubMeshInfo(m, getFaceElementTransformation(m, be), facesInfo, sets.first.second);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing1D(Mesh& m)
{
	auto be_counter{ 0 };
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn)) {
			be_counter++;
		}
	}
	switch (be_counter) {
	case 1:
		assignIndividualTFSFAttsOnePoint1D(m);
		break;
	case 2:
		assignIndividualTFSFAttsTwoPoints1D(m);
		break;
	default:
		throw std::exception("Only one or two TFSF points can be declared in a 1D Mesh.");
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing(Mesh& m)
{
	auto facemap{ m.GetFaceToBdrElMap() };
	for (int be = 0; be < m.GetNBE(); be++) {
		if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldIn) ||
			m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {

			auto be_trans{ getFaceElementTransformation(m, be) };
			auto face_oris{ calculateBaryNormalProduct(m, *be_trans, be) };

			if (m.GetBdrAttribute(be) == static_cast<int>(BdrCond::TotalFieldOut)) {
				face_oris.first  *= -1.0; 
				face_oris.second *= -1.0;
			}

			Array<int> be_vert, el1_face, el1_ori, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			switch (m.Dimension()) {
			case 1:
				m.GetElementVertices(be_trans->Elem1No, el1_face);
				break;
			case 2:
				m.GetElementEdges(be_trans->Elem1No, el1_face, el1_ori);
				break;
			case 3:
				m.GetElementFaces(be_trans->Elem1No, el1_face, el1_ori);
				break;
			default:
				throw std::exception("Incorrect Dimension for mesh in TFSF Child Attribute setting.");
			}
			std::pair<FaceId, IsTF> set_v1;
			for (int f = 0; f < el1_face.Size(); f++) {
				auto fi{ m.GetFaceInformation(f) };
				m.GetFaceVertices(el1_face[f], face_vert);
				face_vert.Sort();
				if (face_vert == be_vert) {
					face_oris.first >= 0.0 ? set_v1 = std::make_pair(f, true) : set_v1 = std::make_pair(f, false);
					break;
				}
			}

			std::pair<FaceId, IsTF> set_v2;
			if (be_trans->Elem2No != NotFound) {
				switch (m.Dimension()) {
				case 1:
					m.GetElementVertices(be_trans->Elem2No, el2_face);
					break;
				case 2:
					m.GetElementEdges(be_trans->Elem2No, el2_face, el2_ori);
					break;
				case 3:
					m.GetElementFaces(be_trans->Elem2No, el2_face, el2_ori);
					break;
				default:
					throw std::exception("Incorrect Dimension for mesh in TFSF Child Attribute setting.");
				}
				for (int f = 0; f < el2_face.Size(); f++) {
					auto fi{ m.GetFaceInformation(f) };
					auto ir = Geometries.GetVertices(Geometry::Type::SQUARE);
					m.GetFaceVertices(el2_face[f], face_vert);
					auto el_faces = m.GetElement(be_trans->Elem2No)->GetFaceVertices(f);
					face_vert.Sort();
					if (face_vert == be_vert) {
						face_oris.second >= 0.0 ? set_v2 = std::make_pair(f, true) : set_v2 = std::make_pair(f, false);
						break;
					}
				}
			}
			else {
				auto set_v2{ std::make_pair(NotFound, false) };
			}
			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.
			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			std::pair<IsTF, IsTF> dirInfo = std::make_pair(set_v1.second, set_v2.second);
			prepareSubMeshInfo(m, be_trans, facesInfo, set_v1.second);
		}
	}
}


MaxwellTransferMap::MaxwellTransferMap(const GridFunction& src,
	const GridFunction& dst)
{
	const FiniteElementSpace* parentfes = nullptr, * subfes = nullptr;
	SubMesh* src_sm = static_cast<SubMesh*>(src.FESpace()->GetMesh());
	subfes = src.FESpace();
	parentfes = dst.FESpace();
	SubMeshUtils::BuildVdofToVdofMap(*subfes,
		*parentfes,
		src_sm->GetFrom(),
		src_sm->GetParentElementIDMap(),
		sub_to_parent_map_);
}

void MaxwellTransferMap::TransferAdd(const GridFunction& src, GridFunction& dst) const
{
	for (int i = 0; i < sub_to_parent_map_.Size(); i++)
	{
		dst(sub_to_parent_map_[i]) += src(i);
	}
}

void MaxwellTransferMap::TransferSub(const GridFunction& src, GridFunction& dst) const
{
	for (int i = 0; i < sub_to_parent_map_.Size(); i++)
	{
		dst(sub_to_parent_map_[i]) -= src(i);
	}
}

};