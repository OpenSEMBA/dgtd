#include "SubMesher.h"

namespace maxwell {

using namespace mfem;
static const int NotFound{ -1 };

FaceElementTransformations* getFaceElementTransformation(Mesh&m, int be) 
{
	switch (m.FaceIsInterior(m.GetBdrElementFaceIndex(be))) {
	case true:
		return m.GetInternalBdrFaceTransformations(be);
	default:
		return m.GetBdrFaceTransformations(be);
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
		throw std::runtime_error("Wrong dimension for TFSF barycenter calculation.");
	}
	res = 0.0;
	for (int v = 0; v < elem_vert.Size(); v++) {
		Vector vertexPos(m.GetVertex(elem_vert[v]), m.Dimension());
		res += vertexPos;
	}
	res /= elem_vert.Size();
	return res;
}

Vector getBarycenterOfFaceElement(Mesh& m, int f)
{
	Array<int> f_elem_vert;
	m.GetFaceVertices(f, f_elem_vert);
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
		throw std::runtime_error("Wrong dimension for TFSF barycenter calculation.");
	}
	res = 0.0;
	for (int v = 0; v < f_elem_vert.Size(); v++) {
		Vector vertexPos(m.GetVertex(f_elem_vert[v]), m.Dimension());
		res += vertexPos;
	}
	res /= f_elem_vert.Size();
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
		throw std::runtime_error("Wrong dimension for TFSF barycenter calculation.");
	}
	for (int i = 0; i < res.Size(); i++) {
		res[i] = b_v[i] - bdr_v[i];
	}
	return res;
}

Vector getNormal(FaceElementTransformations& fet)
{
	Vector res;
	switch (fet.Elem1->GetDimension()) {
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
		throw std::runtime_error("Wrong Dimension for Element in Normal for TFSF orientations.");
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

double calculateCrossBaryVertexSign(Mesh& m, FaceElementTransformations& fet, int be)
{
	auto coord_v0{ m.GetVertex(m.GetBdrElement(be)->GetVertices()[0]) };
	auto coord_v1{ m.GetVertex(m.GetBdrElement(be)->GetVertices()[1]) };

	auto bary_e1{ getBarycenterOfElement(m, fet.Elem1No) };
	auto bary_e2{ getBarycenterOfElement(m, fet.Elem2No) };

	Vector bary_3D(3), vertex_3D(3);
	bary_3D = 0.0; vertex_3D = 0.0;
	for (int i = 0; i < m.Dimension(); i++) {
		bary_3D[i] = bary_e2[i] - bary_e1[i];
		vertex_3D[i] = coord_v1[i] - coord_v0[i];
	}

	auto cross{ crossProduct(bary_3D,vertex_3D) };

	return cross[2];

}

const Vector buildTangent2D(Mesh& m, int be)
{
	auto be_trans{ m.GetBdrElementTransformation(be) };
	auto f{ m.GetFace(m.GetBdrElementFaceIndex(be)) };
	auto v0{ m.GetVertex(f->GetVertices()[0]) };
	auto v1{ m.GetVertex(f->GetVertices()[1]) };
	Vector res(2);
	for (int i = 0; i < res.Size(); i++) {
		res[i] = v1[i] - v0[i];
	}
	return res;
}

const Vector buildNormal3D(Mesh& m, int be)
{
	auto be_trans{ m.GetBdrElementTransformation(be) };
	auto f{ m.GetFace(m.GetBdrElementFaceIndex(be)) };
	auto v0{ m.GetVertex(f->GetVertices()[0]) };
	auto v1{ m.GetVertex(f->GetVertices()[1]) };
	double* v2;
	f->GetGeometryType() == Geometry::SQUARE ? v2 = m.GetVertex(f->GetVertices()[3]) : v2 = m.GetVertex(f->GetVertices()[2]);
	Vector aux_vec1(3), aux_vec2(3);
	for (int i = 0; i < aux_vec1.Size(); i++) {
		aux_vec1[i] = v1[i] - v0[i];
		aux_vec2[i] = v2[i] - v0[i];
	}
	auto res{ crossProduct(aux_vec1, aux_vec2) };
	return res;
}

const std::pair<Vector, Vector> calculateBarycenters(Mesh& m, int be)
{
	auto fe_trans{ getFaceElementTransformation(m, be) };

	if (fe_trans->Elem2No != NotFound) {
		return std::make_pair(getBarycenterOfElement(m, fe_trans->Elem1No), getBarycenterOfElement(m, fe_trans->Elem2No));
	}
	else {
		return std::make_pair(getBarycenterOfElement(m, fe_trans->Elem1No), getBarycenterOfElement(m, m.GetBdrElementFaceIndex(be)));
	}
}

const Vector buildBarycenterPosition(Mesh& m, int be)
{
	auto barys{ calculateBarycenters(m, be) };

	Vector res(barys.first.Size());
	for (auto i{ 0 }; i < res.Size(); ++i) {
		res[i] = barys.second[i] - barys.first[i];
	}
	return res;
}

void restoreElementAttributes(Mesh& m) //Temporary method that has to be reworked when materials are implemented.
{
	for (int e = 0; e < m.GetNE(); e++)
	{
		m.SetAttribute(e, 1);
	}
}


void cleanInvalidSubMeshEntries(std::vector<El2Face>& v)
{
	auto end_tf = std::remove_if(v.begin(), v.end(), [](const auto& i) {
		return i.first == -1;
		});
	v.erase(end_tf, v.end());
}


void setBoundaryAttributesInChild(const Mesh& parent, SubMesh& child, const std::pair<Array<int>, BdrCond>& parent_info)
{
	//if (child.Dimension() == 1) {
	//	for (int e = 0; e < child.GetNE(); e++) {
	//		Array<int> verts(2);
	//		child.GetElementVertices(e, verts);
	//		child.AddBdrPoint(verts[0]);
	//		child.AddBdrPoint(verts[1]);
	//	}
	//}
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int be = 0; be < parent.GetNBE(); be++) {
		if (parent_info.first[parent.GetBdrAttribute(be) - 1] == 1) {
			child.SetBdrAttribute(child_f2bdr_map[map.Find(parent_f2bdr_map.Find(be))], static_cast<int>(parent_info.second));
		}
	}
}

Array<int> getMarkerForSubMesh(const BdrCond& bdrCond, bool isTF)
{
	Array<int> res(1);
	switch (bdrCond) {
	case BdrCond::TotalFieldIn:
		if (isTF) {
			res[0] = SubMeshingMarkers::TotalFieldMarker;
		}
		else {
			res[0] = SubMeshingMarkers::ScatteredFieldMarker;
		}
		break;
	case BdrCond::NearToFarField:
		res[0] = SubMeshingMarkers::NearToFarFieldMarker;
		break;
	}
	return res;
}

SubMesh createSubMeshFromParent(const Mesh& parent, const std::pair<Array<int>, BdrCond>& parent_info, bool isTF = false)
{
	Array<int> marker{ getMarkerForSubMesh(parent_info.second, isTF) };

	auto res{ SubMesh::CreateFromDomain(parent, marker) };
	setBoundaryAttributesInChild(parent, res, parent_info);

	restoreElementAttributes(res);
	res.FinalizeMesh();
	return res;
}

TotalFieldScatteredFieldSubMesher::TotalFieldScatteredFieldSubMesher(const Mesh& m, const Array<int>& marker)
{
	Mesh parent_for_global(m);
	Mesh parent_for_individual(m);

	setGlobalTFSFAttributesForSubMeshing(parent_for_global, marker);

	switch (m.Dimension()) {
	case 1:
		setIndividualTFSFAttributesForSubMeshing1D(parent_for_individual, marker);
		break;
	case 2:
		setIndividualTFSFAttributesForSubMeshing2D(parent_for_individual, marker);
		break;
	default:
		setIndividualTFSFAttributesForSubMeshing3D(parent_for_individual, marker);
		break;
	}

	Array<int> global_att(1); global_att[0] = SubMeshingMarkers::GlobalSubMeshMarker;
	auto global_sm{ SubMesh::CreateFromDomain(parent_for_global, global_att) };
	restoreElementAttributes(global_sm);
	global_sm.FinalizeMesh();
	global_submesh_ = std::make_unique<SubMesh>(global_sm);

	cleanInvalidSubMeshEntries(elem_to_face_tf_);
	cleanInvalidSubMeshEntries(elem_to_face_sf_);

	if (!elem_to_face_tf_.empty()) {
		tf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, std::make_pair(marker, BdrCond::TotalFieldIn), true));
	}

	if (!elem_to_face_sf_.empty()) {
		sf_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(parent_for_individual, std::make_pair(marker, BdrCond::TotalFieldIn), false));
	}

};



void TotalFieldScatteredFieldSubMesher::setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el1_is_tf)
{
	if (trans->Elem2No != NotFound) {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalFieldMarker);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::ScatteredFieldMarker);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredFieldMarker);
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::TotalFieldMarker);
		}
	}
	else {
		if (el1_is_tf) {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::TotalFieldMarker);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::ScatteredFieldMarker);
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


void TotalFieldScatteredFieldSubMesher::setGlobalTFSFAttributesForSubMeshing(Mesh& m, const Array<int>& marker)
{

	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {
			auto be_trans{ getFaceElementTransformation(m, be)};
			if (be_trans->Elem2No != NotFound) {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::GlobalSubMeshMarker);
				m.GetElement(be_trans->Elem2No)->SetAttribute(SubMeshingMarkers::GlobalSubMeshMarker);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
				elems_for_global_submesh_.push_back(be_trans->Elem2No);
			}
			else {
				m.GetElement(be_trans->Elem1No)->SetAttribute(SubMeshingMarkers::GlobalSubMeshMarker);
				elems_for_global_submesh_.push_back(be_trans->Elem1No);
			}
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

void TotalFieldScatteredFieldSubMesher::assignIndividualTFSFAttsOnePoint1D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {
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

void TotalFieldScatteredFieldSubMesher::assignIndividualTFSFAttsTwoPoints1D(Mesh& m, const Array<int>& marker)
{
	auto flag{ false };
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {
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

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing1D(Mesh& m, const Array<int>& marker)
{
	auto be_counter{ 0 };
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {
			be_counter++;
		}
	}
	switch (be_counter) {
	case 1:
		assignIndividualTFSFAttsOnePoint1D(m, marker);
		break;
	case 2:
		assignIndividualTFSFAttsTwoPoints1D(m, marker);
		break;
	default:
		throw std::runtime_error("Only one or two TFSF points can be declared in a 1D Mesh.");
	}
}

double buildCrossCoefficient(const Vector& bary_vec, const Vector& tang_vec)
{
	Vector cross_first(3), cross_sec(3);
	cross_first[0] = bary_vec[0];
	cross_first[1] = bary_vec[1];
	cross_first[2] = 0.0;
	cross_sec[0] = tang_vec[0];
	cross_sec[1] = tang_vec[1];
	cross_sec[2] = 0.0;
	auto cross = crossProduct(cross_first, cross_sec);
	return cross[2];
}

double buildFaceOrientation(Mesh& mesh, int be)
{
	switch (mesh.Dimension())
	{
	case 2:
		return buildCrossCoefficient(
			buildBarycenterPosition(mesh, be),
			buildTangent2D(mesh, be));
	case 3:
		return mfem::InnerProduct(
			buildBarycenterPosition(mesh, be),
			buildNormal3D(mesh, be));
	default:
		throw std::runtime_error("Method only supports 2D and 3D.");
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing2D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			auto fe_trans{ getFaceElementTransformation(m, be) };

			Array<int> be_vert, el1_face, el1_ori, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			m.GetElementEdges(fe_trans->Elem1No, el1_face, el1_ori);

			std::pair<FaceId, IsTF> set_v1;
			for (int f = 0; f < el1_face.Size(); f++) {
				auto fi{ m.GetFaceInformation(f) };
				m.GetEdgeVertices(el1_face[f], face_vert);
				face_vert.Sort();
				if (face_vert == be_vert) {
					face_ori >= 0.0 ? set_v1 = std::make_pair(f, false) : set_v1 = std::make_pair(f, true);
					break;
				}
			}

			std::pair<FaceId, IsTF> set_v2;
			if (fe_trans->Elem2No != NotFound) {
				m.GetElementEdges(fe_trans->Elem2No, el2_face, el2_ori);
				for (int f = 0; f < el2_face.Size(); f++) {
					auto fi{ m.GetFaceInformation(f) };
					auto ir = Geometries.GetVertices(Geometry::Type::SQUARE);
					m.GetFaceVertices(el2_face[f], face_vert);
					auto el_faces = m.GetElement(fe_trans->Elem2No)->GetEdgeVertices(f);
					face_vert.Sort();
					if (face_vert == be_vert) {
						face_ori >= 0.0 ? set_v2 = std::make_pair(f, true) : set_v2 = std::make_pair(f, false);
						break;
					}
				}
			}
			else {
				auto set_v2{ std::make_pair(NotFound, false) };
			}
			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.
			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			prepareSubMeshInfo(m, fe_trans, facesInfo, set_v1.second);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing3D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			Array<int> be_vert, el1_face, el1_ori, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m,be) };
			m.GetElementFaces(fe_trans->Elem1No, el1_face, el1_ori);

			std::pair<FaceId, IsTF> set_v1;
			for (int f = 0; f < el1_face.Size(); f++) {
				auto fi{ m.GetFaceInformation(f) };
				m.GetFaceVertices(el1_face[f], face_vert);
				face_vert.Sort();
				if (face_vert == be_vert) {
					face_ori >= 0.0 ? set_v1 = std::make_pair(f, false) : set_v1 = std::make_pair(f, true);
					break;
				}
			}

			std::pair<FaceId, IsTF> set_v2;
			if (fe_trans->Elem2No != NotFound) {

				m.GetElementFaces(fe_trans->Elem2No, el2_face, el2_ori);

				for (int f = 0; f < el2_face.Size(); f++) {
					auto fi{ m.GetFaceInformation(f) };
					auto ir = Geometries.GetVertices(Geometry::Type::SQUARE);
					m.GetFaceVertices(el2_face[f], face_vert);
					auto el_faces = m.GetElement(fe_trans->Elem2No)->GetFaceVertices(f);
					face_vert.Sort();
					if (face_vert == be_vert) {
						face_ori >= 0.0 ? set_v2 = std::make_pair(f, true) : set_v2 = std::make_pair(f, false);
						break;
					}
				}
			}
			else {
				auto set_v2{ std::make_pair(NotFound, false) };
			}
			//be_vert is counterclockwise, that is our convention to designate which element will be TF. The other element will be SF.
			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			prepareSubMeshInfo(m, fe_trans, facesInfo, set_v1.second);
		}
	}

}

NearToFarFieldSubMesher::NearToFarFieldSubMesher(const Mesh& m, const FiniteElementSpace& fes, const Array<int>& marker)
{

	original_ = std::make_unique<Mesh>(m);

	switch (original_->SpaceDimension()) {
	case 2:
		setIndividualNTFFAttributesForSubMeshing2D(*original_.get(), marker);
		break;
	case 3:
		setIndividualNTFFAttributesForSubMeshing3D(*original_.get(), marker);
		break;
	default:
		throw std::runtime_error("NearToFarField can only be applied to 2D or 3D meshes.");
	}

	if (!elem_to_face_ntff_.empty()) {
		ntff_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(*original_.get(), std::make_pair(marker, BdrCond::NearToFarField)));
	}
	else {
		throw std::runtime_error("ntff submesh is empty, check orientation of curves/faces.");
	}
}

void NearToFarFieldSubMesher::setIndividualNTFFAttributesForSubMeshing2D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			Array<int> be_vert, el_face, el_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m, be) };

			std::pair<FaceId, IsTF> el_info;
			if (fe_trans->Elem2No != NotFound) {

				m.GetElementEdges(fe_trans->Elem1No, el_face, el_ori);

				for (int f = 0; f < el_face.Size(); f++) {
					auto fi{ m.GetFaceInformation(f) };
					m.GetFaceVertices(el_face[f], face_vert);
					face_vert.Sort();
					if (face_vert == be_vert) {
						face_ori >= 0.0 ? el_info = std::make_pair(f, false) : el_info = std::make_pair(f, true);
						break;
					}
				}
			}
			else {
				std::string error{ "Element 2 has not been found for boundary element " + std::to_string(be) + ", verify NtFF orientations on the mesh adapt to the convention." };
				throw std::runtime_error(error.c_str());
			}
			//Our convention is based on the inner product between a vector that joins the barycenters of the elements (going from elem1 to elem2)
			//and the normal vector on the face, if it's positive, we designate it as TF. The other element will be SF.
			prepareSubMeshInfo(m, fe_trans, el_info.first, el_info.second);
		}
	}

}

void NearToFarFieldSubMesher::setIndividualNTFFAttributesForSubMeshing3D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			Array<int> be_vert, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m, be) };

			std::pair<FaceId, IsTF> el2_info;
			if (fe_trans->Elem2No != NotFound) {

				m.GetElementFaces(fe_trans->Elem2No, el2_face, el2_ori);

				for (int f = 0; f < el2_face.Size(); f++) {
					auto fi{ m.GetFaceInformation(f) };
					m.GetFaceVertices(el2_face[f], face_vert);
					face_vert.Sort();
					if (face_vert == be_vert) {
						face_ori >= 0.0 ? el2_info = std::make_pair(f, false) : el2_info = std::make_pair(f, true);
						break;
					}
				}
			}
			else {
				std::string error{ "Element 2 has not been found for boundary element " + std::to_string(be) + ", verify NtFF orientations on the mesh adapt to the convention." };
				throw std::runtime_error(error.c_str());
			}
			//Our convention is based on the inner product between a vector that joins the barycenters of the elements (going from elem1 to elem2)
			//and the normal vector on the face, if it's positive, we designate it as TF. The other element will be SF.
			prepareSubMeshInfo(m, fe_trans, el2_info.first, el2_info.second);
		}
	}
}

void NearToFarFieldSubMesher::prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, int faceId, bool el_is_ntff)
{
	setAttributeForTagging(m, trans, el_is_ntff);
	storeElementToFaceInformation(trans, faceId, el_is_ntff);
}

void NearToFarFieldSubMesher::setAttributeForTagging(Mesh& m, const FaceElementTransformations* trans, bool el_is_ntff)
{
	if (trans->Elem2No != NotFound) {
		if (el_is_ntff) {
			m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::NearToFarFieldMarker);
		}
		else {
			m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::NearToFarFieldMarker);
		}
	}
}

void NearToFarFieldSubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, int faceId, bool el_is_ntff)
{
	if (faceId != NotFound) {
		if (el_is_ntff) {
			elem_to_face_ntff_.push_back(std::make_pair(trans->Elem2No, faceId));
		}
		else {
			elem_to_face_ntff_.push_back(std::make_pair(trans->Elem1No, faceId));
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
