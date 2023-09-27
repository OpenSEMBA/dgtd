#pragma once

#include <mfem.hpp>

using FaceId = int;
using ElementId = int;
using Attribute = int;
using IsCCW = bool;
using El2Face = std::pair<ElementId, FaceId>;
using Face2Dir = std::pair<FaceId, IsCCW>;

namespace maxwell {

using namespace mfem;

class TotalFieldScatteredFieldSubMesher
{
public:

	TotalFieldScatteredFieldSubMesher(){};
	TotalFieldScatteredFieldSubMesher(const Mesh&);

	SubMesh* getTFSubMesh() { return tf_mesh_.get(); }
	SubMesh* getSFSubMesh() { return sf_mesh_.get(); }
	const SubMesh* getTFConstSubMesh() const { return tf_mesh_.get(); }
	const SubMesh* getSFConstSubMesh() const { return sf_mesh_.get(); }

private:

	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void setBoundaryAttributesInChild(const Mesh& parent, SubMesh& child);
	void storeElementToFaceInformation(const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf);
	void prepareSubMeshInfo(Mesh& m, const FaceElementTransformations* trans, const std::pair<int, int> facesId, bool el1_is_tf);
	void setTFSFAttributesForSubMeshing(Mesh&);
	void restoreElementAttributes(Mesh& m);

	std::pair<FaceId, IsCCW> getFaceAndDirOnVertexIteration(const Element*, const Array<int>& verts, const Array<int>& be_verts);

	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;

	std::unique_ptr<SubMesh> tf_mesh_;
	std::unique_ptr<SubMesh> sf_mesh_;

};

class MaxwellTransferMap
{
public:
	
	MaxwellTransferMap(const GridFunction& src, const GridFunction& dst);
	
	void TransferAdd(const GridFunction& src, GridFunction& dst) const;

private:

	Array<int> sub_to_parent_map_;
	mutable Vector z_;

};

};