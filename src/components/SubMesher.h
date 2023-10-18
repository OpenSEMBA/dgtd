#pragma once

#include <mfem.hpp>
#include "math/Calculus.h"
#include "Types.h"

using FaceId = int;
using ElementId = int;
using Attribute = int;
using IsTF = bool;
using El2Face = std::pair<ElementId, FaceId>;
using Face2Dir = std::pair<FaceId, IsTF>;
using SetPairs = std::pair<std::pair<FaceId, IsTF>, std::pair<FaceId, IsTF>>;

namespace maxwell {

using namespace mfem;

class TotalFieldScatteredFieldSubMesher
{
public:

	TotalFieldScatteredFieldSubMesher(){};
	TotalFieldScatteredFieldSubMesher(const Mesh&);

	SubMesh* getTFSubMesh() { return tf_mesh_.get(); }
	SubMesh* getSFSubMesh() { return sf_mesh_.get(); }
	SubMesh* getGlobalTFSFSubMesh() { return global_submesh_.get(); }
	const SubMesh* getTFConstSubMesh() const { return tf_mesh_.get(); }
	const SubMesh* getSFConstSubMesh() const { return sf_mesh_.get(); }
	const SubMesh* getGlobalTFSFConstSubMesh() const { return global_submesh_.get(); }

private:

	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void setBoundaryAttributesInChild1D(const Mesh& parent, SubMesh& child);
	void setBoundaryAttributesInChild2D(const Mesh& parent, SubMesh& child);
	void setBoundaryAttributesInChild3D(const Mesh& parent, SubMesh& child);
	void storeElementToFaceInformation(const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	void prepareSubMeshInfo(Mesh& m,   const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	void setGlobalTFSFAttributesForSubMeshing(Mesh&);
	
	SetPairs twoPointAssignator(Mesh&, int be, bool flag);
	void assignIndividualTFSFAttsOnePoint1D(Mesh&);
	void assignIndividualTFSFAttsTwoPoints1D(Mesh&);
	void setIndividualTFSFAttributesForSubMeshing1D(Mesh&);
	
	void setIndividualTFSFAttributesForSubMeshing(Mesh&);
	void setIndividualTFSFAttributesForSubMeshing2D(Mesh&);

	void setIndividualTFSFAttributesForSubMeshing3D(Mesh&);

	void restoreElementAttributes(Mesh&);
	FaceElementTransformations* getFaceElementTransformation(Mesh&, int bdr_el_no);
	SubMesh TotalFieldScatteredFieldSubMesher::createSubMeshFromParent(const Mesh&, bool isTF);

	Face2Dir getFaceAndDirOnVertexIteration2D(const Element*, const Array<int>& verts, const Array<int>& be_verts);
	Face2Dir getFaceAndDirOnVertexIteration3D(Mesh& m, int be);

	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;
	std::vector<ElementId> elems_for_global_submesh_;

	std::unique_ptr<SubMesh> tf_mesh_;
	std::unique_ptr<SubMesh> sf_mesh_;
	std::unique_ptr<SubMesh> global_submesh_;

};

class MaxwellTransferMap
{
public:
	
	MaxwellTransferMap(const GridFunction& src, const GridFunction& dst);
	
	void TransferAdd(const GridFunction& src, GridFunction& dst) const;
	void TransferSub(const GridFunction& src, GridFunction& dst) const;

private:

	Array<int> sub_to_parent_map_;
	mutable Vector z_;

};

};