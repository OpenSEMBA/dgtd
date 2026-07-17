#include "SubMesher.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>

namespace maxwell {

using namespace mfem;
static const int NotFound{ -1 };

static void setVectorSizeForDim(int dim, Vector& res)
{
	switch (dim) {
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
}

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const ParFiniteElementSpace& fes)
{
	Array<int> res(fes.GetMesh()->bdr_attributes.Max());
	res = 0;
	for (const int t : tags) {
		res[t - 1] = 1;
	}
	return res;
}

FaceElementTransformations* getFaceElementTransformation(Mesh&m, int be)
{
	if (auto* tr = m.GetInternalBdrFaceTransformations(be)) {
		return tr;
	}
	return m.GetBdrFaceTransformations(be);
}

Vector getBarycenterOfElement(Mesh& m, int e)
{
	auto elem{ m.GetElement(e) };
	Array<int> elem_vert(elem->GetNVertices());
	elem->GetVertices(elem_vert);
	Vector res;
	setVectorSizeForDim(m.Dimension(), res);
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
	setVectorSizeForDim(m.Dimension(), res);
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
	setVectorSizeForDim(b_v.Size(), res);
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

Vector getNormalAtRefPoint(FaceElementTransformations& fet, const IntegrationPoint& face_ip)
{
	Vector res;
	int dim = fet.Elem1->GetDimension();
	switch (dim) {
	case 1:
		res.SetSize(1);
		res(0) = 1.0;
		break;
	case 2:
	case 3:
		res.SetSize(dim);
		fet.SetIntPoint(&face_ip);
		CalcOrtho(fet.Jacobian(), res);
		break;
	default:
		throw std::runtime_error("Wrong dimension in getNormalAtRefPoint.");
	}
	return res;
}

std::vector<Vector> computePerDofFaceNormals(
	FaceElementTransformations& fet,
	const FiniteElementSpace& fes,
	const std::vector<int>& elem1FaceDofLocalIndices)
{
	int dim = fet.Elem1->GetDimension();
	int nDofs = static_cast<int>(elem1FaceDofLocalIndices.size());
	std::vector<Vector> normals(nDofs);

	if (dim <= 1) {
		for (int i = 0; i < nDofs; i++) {
			normals[i].SetSize(1);
			normals[i](0) = 1.0;
		}
		return normals;
	}

	// Get Elem1's finite element and its DOF reference-space positions.
	const FiniteElement* fe = fes.GetFE(fet.Elem1No);
	const IntegrationRule& elemNodes = fe->GetNodes();

	// Get the face's geometry and build a fine integration rule to cover
	// all DOF locations on the face. Use the face FE's nodes if available,
	// otherwise use a high-order rule.
	int faceGeomType = fet.GetGeometryType();
	int order = fe->GetOrder();
	const IntegrationRule& faceIR =
		IntRules.Get(faceGeomType, 2 * order);

	// For each face integration point, transform to Elem1-ref-space via Loc1
	// and compute the physical position. Build a lookup of face ref points
	// mapped to Elem1 reference coords.
	// We'll match element DOFs to face points by proximity in element ref space.

	for (int di = 0; di < nDofs; di++) {
		int localDof = elem1FaceDofLocalIndices[di];
		const IntegrationPoint& dofIP = elemNodes.IntPoint(localDof);

		// Find the face reference point that maps closest to this DOF's
		// element reference position, via Loc1.
		double bestDist = 1e30;
		int bestFI = 0;
		for (int fi = 0; fi < faceIR.GetNPoints(); fi++) {
			IntegrationPoint elemIP;
			fet.Loc1.Transform(faceIR.IntPoint(fi), elemIP);
			double dx = elemIP.x - dofIP.x;
			double dy = elemIP.y - dofIP.y;
			double dz = elemIP.z - dofIP.z;
			double dist = dx*dx + dy*dy + dz*dz;
			if (dist < bestDist) {
				bestDist = dist;
				bestFI = fi;
			}
		}

		normals[di] = getNormalAtRefPoint(fet, faceIR.IntPoint(bestFI));
	}

	return normals;
}

std::pair<double, double> calculateBaryNormalProduct(ParMesh& m, FaceElementTransformations& fet, int be)
{
	auto bdr_vertices{ m.GetBdrElement(be)->GetVertices() };
	auto bdr_vert{ m.GetVertex(bdr_vertices[0]) };
	auto normal{ getNormal(fet) };

	auto bary1{ getBarycenterOfElement(m, fet.Elem1No) };
	auto v1{ subtract(bdr_vert,bary1) };
	auto d1{ mfem::InnerProduct(v1, normal) };

	auto d2 = 0.0;
	if (fet.Elem2No >= 0) {
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

	auto baris = calculateBarycenters(m, be);

	Vector bary_3D(3), vertex_3D(3);
	bary_3D = 0.0; vertex_3D = 0.0;
	for (int i = 0; i < m.Dimension(); i++) {
		bary_3D[i] = baris.second[i] - baris.first[i];
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

	if (fe_trans->Elem2No >= 0) {
		return std::make_pair(getBarycenterOfElement(m, fe_trans->Elem1No), getBarycenterOfElement(m, fe_trans->Elem2No));
	}
	else {
		return std::make_pair(getBarycenterOfElement(m, fe_trans->Elem1No), getBarycenterOfFaceElement(m, m.GetBdrElementFaceIndex(be)));
	}
}

const Vector buildElem1ToElem2BarycenterVector(Mesh& m, int be)
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
	auto parent_f2bdr_map{ parent.GetFaceToBdrElMap() };
	auto child_f2bdr_map{ child.GetFaceToBdrElMap() };
	auto map{ SubMeshUtils::BuildFaceMap(parent, child, child.GetParentElementIDMap()) };
	for (int be = 0; be < parent.GetNBE(); be++) {
		if (parent_info.first[parent.GetBdrAttribute(be) - 1] == 1) {
			const int parent_face = parent_f2bdr_map.Find(be);
			if (parent_face < 0) continue;
			const int child_face = map.Find(parent_face);
			if (child_face < 0) continue; // face not in child SubMesh (e.g. partition boundary)
			const int child_bdr = child_f2bdr_map[child_face];
			if (child_bdr < 0) continue;
			child.SetBdrAttribute(child_bdr, static_cast<int>(parent_info.second));
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
	if (trans->Elem2No >= 0) {
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
			if (be_trans->Elem2No >= 0) {
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
		if (be_trans->Elem1No == 0 && be_trans->Elem2No < 0) {
			set_e1 = std::make_pair(0, true);
			set_e2 = std::make_pair(NotFound, false);
		}
		else {
			set_e1 = std::make_pair(1, false);
			set_e2 = std::make_pair(0, true);
		}
		break;
	case true: //If second boundary we find
		if (be_trans->Elem1No == m.GetNE() - 1 && be_trans->Elem2No < 0) {
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
			buildElem1ToElem2BarycenterVector(mesh, be),
			buildTangent2D(mesh, be));
	case 3:
		return mfem::InnerProduct(
			buildElem1ToElem2BarycenterVector(mesh, be),
			buildNormal3D(mesh, be));
	default:
		throw std::runtime_error("Method only supports 2D and 3D.");
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing2D(Mesh& m, const Array<int>& marker)
{

	/*
							 @   MFEM numbering convention is read from
							@@@@   the GMSH element order, thus we need to ensure
						  @@ @  @@   when designing the mesh that the numbering order
						@@   @    @   of the elements and the 2D vector along the face align
					  @@     @     @@   with out TFSF designation intent
					 @       @       @@
				   @@        @         @@  
				 @@          A           @@  
				@            |            @@
			  @@    SF       |      TF      @@
			@@               |                @
		   @@    elem 1  ----+---> elem 2      @
			 @@   i.e. 75    |    i.e. 90   @@
			  @@             |             @@
				@@           |           @@
				  @@         @          @
					@@       @        @@
					  @@     @      @@
					   @@    @     @
						 @@  @   @@
						   @ @ @@
							@@@
							 @
	 */

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
			if (fe_trans->Elem2No >= 0) {
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

			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			prepareSubMeshInfo(m, fe_trans, facesInfo, set_v1.second);
		}
	}
}

void TotalFieldScatteredFieldSubMesher::setIndividualTFSFAttributesForSubMeshing3D(Mesh& m, const Array<int>& marker)
{
	// Compute centroid of all TFSF boundary vertices.
	// For a convex TFSF surface the centroid lies inside the TF region,
	// so the element whose barycenter is closer to it is the TF element.
	// This is invariant to face-normal orientation and handles corner
	// elements (adjacent to multiple TFSF faces) consistently.
	Vector tfsf_center(3);
	tfsf_center = 0.0;
	int n_verts = 0;
	std::set<int> counted;
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] != 1) continue;
		Array<int> verts;
		m.GetBdrElementVertices(be, verts);
		for (int i = 0; i < verts.Size(); i++) {
			if (counted.insert(verts[i]).second) {
				auto* v = m.GetVertex(verts[i]);
				for (int d = 0; d < 3; d++) tfsf_center[d] += v[d];
				n_verts++;
			}
		}
	}
	if (n_verts > 0) tfsf_center /= double(n_verts);

	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			Array<int> be_vert, el1_face, el1_ori, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m,be) };
			m.GetElementFaces(fe_trans->Elem1No, el1_face, el1_ori);

			// Determine TF/SF by centroid distance: closer to center = TF.
			Vector bary1 = getBarycenterOfElement(m, fe_trans->Elem1No);
			double dist1_sq = 0.0;
			for (int d = 0; d < 3; d++) {
				double diff = bary1[d] - tfsf_center[d];
				dist1_sq += diff * diff;
			}

			bool elem1_is_tf;
			if (fe_trans->Elem2No >= 0) {
				Vector bary2 = getBarycenterOfElement(m, fe_trans->Elem2No);
				double dist2_sq = 0.0;
				for (int d = 0; d < 3; d++) {
					double diff = bary2[d] - tfsf_center[d];
					dist2_sq += diff * diff;
				}
				elem1_is_tf = dist1_sq < dist2_sq;
			}
			else {
				// Boundary face with no Elem2: compare element to face barycenter
				Vector face_bary = getBarycenterOfFaceElement(m, m.GetBdrElementFaceIndex(be));
				double face_dist_sq = 0.0;
				for (int d = 0; d < 3; d++) {
					double diff = face_bary[d] - tfsf_center[d];
					face_dist_sq += diff * diff;
				}
				elem1_is_tf = dist1_sq < face_dist_sq;
			}

			std::pair<FaceId, IsTF> set_v1;
			for (int f = 0; f < el1_face.Size(); f++) {
				m.GetFaceVertices(el1_face[f], face_vert);
				face_vert.Sort();
				if (face_vert == be_vert) {
					set_v1 = std::make_pair(f, elem1_is_tf);
					break;
				}
			}

			std::pair<FaceId, IsTF> set_v2;
			if (fe_trans->Elem2No >= 0) {

				m.GetElementFaces(fe_trans->Elem2No, el2_face, el2_ori);

				for (int f = 0; f < el2_face.Size(); f++) {
					m.GetFaceVertices(el2_face[f], face_vert);
					face_vert.Sort();
					if (face_vert == be_vert) {
						set_v2 = std::make_pair(f, !elem1_is_tf);
						break;
					}
				}
			}
			else {
				auto set_v2{ std::make_pair(NotFound, false) };
			}

			std::pair<FaceId, FaceId> facesInfo = std::make_pair(set_v1.first, set_v2.first);
			prepareSubMeshInfo(m, fe_trans, facesInfo, set_v1.second);
		}
	}

}

NearToFarFieldSubMesher::NearToFarFieldSubMesher(const Mesh& m, const ParFiniteElementSpace& fes, const Array<int>& marker)
{
	Mesh tmesh(m);
	original_ = std::make_unique<Mesh>(m);

	switch (original_->SpaceDimension()) {
	case 2:
		setSurfaceAttributesForSubMesh2D(*original_.get(), marker);
		break;
	case 3:
		setSurfaceAttributesForSubMesh3D(*original_.get(), marker);
		break;
	default:
		throw std::runtime_error("NearToFarField can only be applied to 2D or 3D meshes.");
	}

	if (!elem_to_face_ntff_.empty()) {
		ntff_mesh_ = std::make_unique<SubMesh>(createSubMeshFromParent(*original_.get(), std::make_pair(marker, BdrCond::NearToFarField)));
	}
	// Ranks with no local NTF surface leave ntff_mesh_ empty.
}

void NearToFarFieldSubMesher::setSurfaceAttributesForSubMesh2D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			Array<int> be_vert, el_face, el_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m, be) };

			std::pair<FaceId, IsTF> el_info;
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
			//Our convention is based on the inner product between a vector that joins the barycenters of the elements (going from elem1 to elem2)
			//and the normal vector on the face, if it's positive, we designate it as TF. The other element will be SF.
			prepareSubMeshInfo(m, fe_trans, el_info.first, el_info.second);
		}
	}

}

void NearToFarFieldSubMesher::setSurfaceAttributesForSubMesh3D(Mesh& m, const Array<int>& marker)
{
	for (int be = 0; be < m.GetNBE(); be++) {
		if (marker[m.GetBdrAttribute(be) - 1] == 1) {

			auto face_ori{ buildFaceOrientation(m, be) };

			Array<int> be_vert, el2_face, el2_ori, face_vert;
			m.GetBdrElementVertices(be, be_vert);
			be_vert.Sort();

			auto fe_trans{ getFaceElementTransformation(m, be) };

			std::pair<FaceId, IsTF> el2_info;
			if (fe_trans->Elem2No >= 0) {

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
				// Partition boundary: Elem2 is on another rank; use Elem1 to find the face index.
				Array<int> el1_face, el1_ori;
				m.GetElementFaces(fe_trans->Elem1No, el1_face, el1_ori);
				for (int f = 0; f < el1_face.Size(); f++) {
					m.GetFaceVertices(el1_face[f], face_vert);
					face_vert.Sort();
					if (face_vert == be_vert) {
						el2_info = std::make_pair(f, false); // el_is_ntff=false → tag Elem1
						break;
					}
				}
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
	if (el_is_ntff && trans->Elem2No >= 0) {
		m.GetElement(trans->Elem2No)->SetAttribute(SubMeshingMarkers::NearToFarFieldMarker);
	}
	else if (!el_is_ntff) {
		m.GetElement(trans->Elem1No)->SetAttribute(SubMeshingMarkers::NearToFarFieldMarker);
	}
}

void NearToFarFieldSubMesher::storeElementToFaceInformation(const FaceElementTransformations* trans, int faceId, bool el_is_ntff)
{
	if (faceId != NotFound) {
		if (el_is_ntff && trans->Elem2No >= 0) {
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

VolumetricRegionSubMesher::VolumetricRegionSubMesher(
	const Mesh& parent,
	const RegionTagSet& vacuum_tags,
	const RegionTagSet& pml_tags)
{
	Mesh parent_copy(parent);
	buildRegionMarkers(parent_copy, vacuum_tags, pml_tags);

	// Detect interface on an untouched parent copy before any submesh extraction,
	// which may alter parent-side bookkeeping in MFEM internals.
	detectVacuumPMLInterface(parent_copy, vacuum_tags, pml_tags);

	if (vacuum_marker_.Size() != 0 && vacuum_marker_.Sum() > 0) {
		auto sub = SubMesh::CreateFromDomain(parent_copy, vacuum_marker_);
		restoreElementAttributes(sub);
		sub.FinalizeMesh();
		vacuum_mesh_ = std::make_unique<SubMesh>(sub);
	}

	if (pml_marker_.Size() != 0 && pml_marker_.Sum() > 0) {
		auto sub = SubMesh::CreateFromDomain(parent_copy, pml_marker_);
		restoreElementAttributes(sub);
		sub.FinalizeMesh();
		pml_mesh_ = std::make_unique<SubMesh>(sub);
	}
}

void VolumetricRegionSubMesher::buildRegionMarkers(
	const Mesh& parent,
	const RegionTagSet& vacuum_tags,
	const RegionTagSet& pml_tags)
{
	const int max_attr = parent.attributes.Max();
	vacuum_marker_.SetSize(max_attr);
	vacuum_marker_ = 0;
	pml_marker_.SetSize(max_attr);
	pml_marker_ = 0;

	for (const int tag : vacuum_tags) {
		if (tag <= 0 || tag > max_attr) {
			throw std::runtime_error("Vacuum region tag is out of mesh attribute range.");
		}
		if (parent.attributes.Find(tag) == -1) {
			throw std::runtime_error("Vacuum region tag is not present in mesh attributes.");
		}
		vacuum_marker_[tag - 1] = 1;
	}

	for (const int tag : pml_tags) {
		if (tag <= 0 || tag > max_attr) {
			throw std::runtime_error("PML region tag is out of mesh attribute range.");
		}
		if (parent.attributes.Find(tag) == -1) {
			throw std::runtime_error("PML region tag is not present in mesh attributes.");
		}
		pml_marker_[tag - 1] = 1;
	}
}

void VolumetricRegionSubMesher::detectVacuumPMLInterface(
	Mesh& parent,
	const RegionTagSet& vacuum_tags,
	const RegionTagSet& pml_tags)
{
	const int bdr_max = parent.bdr_attributes.Max();
	interface_marker_.SetSize(bdr_max);
	interface_marker_ = 0;
	interface_faces_.clear();

	std::unordered_set<int> unique_faces;
	const auto face2bdr = parent.GetFaceToBdrElMap();
	int faces_checked = 0;
	int two_sided = 0;
	int crossings = 0;
	int vac_vac_pairs = 0;
	int pml_pml_pairs = 0;
	int vac_pml_pairs = 0;
	int other_pairs = 0;

	auto is_cross = [&](int attr1, int attr2) {
		const bool e1_vac = vacuum_tags.find(attr1) != vacuum_tags.end();
		const bool e2_vac = vacuum_tags.find(attr2) != vacuum_tags.end();
		const bool e1_pml = pml_tags.find(attr1) != pml_tags.end();
		const bool e2_pml = pml_tags.find(attr2) != pml_tags.end();
		return (e1_vac && e2_pml) || (e2_vac && e1_pml);
	};

	for (int f = 0; f < parent.GetNumFaces(); ++f) {
		++faces_checked;

		const auto* ft = parent.GetFaceElementTransformations(f);
		if (!ft || ft->Elem1No < 0 || ft->Elem2No < 0) {
			continue;
		}

		++two_sided;

		const int e1_attr = parent.GetAttribute(ft->Elem1No);
		const int e2_attr = parent.GetAttribute(ft->Elem2No);
		const bool e1_vac = vacuum_tags.find(e1_attr) != vacuum_tags.end();
		const bool e2_vac = vacuum_tags.find(e2_attr) != vacuum_tags.end();
		const bool e1_pml = pml_tags.find(e1_attr) != pml_tags.end();
		const bool e2_pml = pml_tags.find(e2_attr) != pml_tags.end();

		if (e1_vac && e2_vac) {
			++vac_vac_pairs;
		} else if (e1_pml && e2_pml) {
			++pml_pml_pairs;
		} else if ((e1_vac && e2_pml) || (e2_vac && e1_pml)) {
			++vac_pml_pairs;
		} else {
			++other_pairs;
		}

		if (!is_cross(e1_attr, e2_attr)) {
			continue;
		}

		++crossings;
		unique_faces.insert(f);

		if (f < face2bdr.Size() && face2bdr[f] != -1) {
			const int bdr_attr = parent.GetBdrAttribute(face2bdr[f]);
			if (bdr_attr > 0 && bdr_attr <= bdr_max) {
				interface_marker_[bdr_attr - 1] = 1;
			}
		}
	}

	interface_faces_.assign(unique_faces.begin(), unique_faces.end());
	std::sort(interface_faces_.begin(), interface_faces_.end());

	if (Mpi::WorldRank() == 0) {
		std::cout << "[VolPML InterfaceDetect] face-sweep: total_faces=" << parent.GetNumFaces()
			<< ", checked=" << faces_checked
			<< ", two_sided=" << two_sided
			<< ", crossings=" << crossings
			<< ", vac_vac_pairs=" << vac_vac_pairs
			<< ", pml_pml_pairs=" << pml_pml_pairs
			<< ", vac_pml_pairs=" << vac_pml_pairs
			<< ", other_pairs=" << other_pairs
			<< ", unique_faces=" << interface_faces_.size()
			<< ", marker_sum=" << interface_marker_.Sum()
			<< std::endl;
	}
}

};
