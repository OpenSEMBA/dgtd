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

private:

	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void storeElementToFaceInformation(const FaceElementTransformations*, const FaceId, bool el1_is_tf);
	void prepareSubMeshInfo(Mesh&, const FaceElementTransformations*, const FaceId, bool el1_is_tf);
	void setTFSFAttributesForSubMeshing(Mesh&);

	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;
	std::vector<El2Face> elem_to_att_tf_;
	std::vector<El2Face> elem_to_att_sf_;

};

};