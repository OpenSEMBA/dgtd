#pragma once

#include <mfem.hpp>

using FaceId = int;
using ElementId = int;
using Attribute = int;
using El2Face = std::pair<ElementId, FaceId>;
using El2Att = std::pair<ElementId, Attribute>;

namespace maxwell {

using namespace mfem;

class TotalFieldScatteredField_SubMesher
{
public:

	TotalFieldScatteredField_SubMesher(const Mesh&);

	Mesh& getTFMesh() { return tf_mesh_; }
	Mesh& getSFMesh() { return sf_mesh_; }
	const Mesh& getTFConstMesh() { return tf_mesh_; }
	const Mesh& getSFConstMesh() { return sf_mesh_; }

private:

	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void setBoundaryAttributesInChild(const Mesh& parent, SubMesh& child);
	void storeElementToFaceInformation(const FaceElementTransformations*, const FaceId, bool el1_is_tf);
	void prepareSubMeshInfo(Mesh&, const FaceElementTransformations*, const FaceId, bool el1_is_tf);
	void setTFSFAttributesForSubMeshing(Mesh&);

	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;
	std::vector<El2Face> elem_to_att_tf_;
	std::vector<El2Face> elem_to_att_sf_;

	Mesh tf_mesh_;
	Mesh sf_mesh_;

};

};