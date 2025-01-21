#pragma once

#include <mfem.hpp>
#include "math/Calculus.h"
#include "Types.h"


namespace maxwell {

using Attribute = int;
using IsTF = bool;
using Face2Dir = std::pair<FaceId, IsTF>;
using SetPairs = std::pair<std::pair<FaceId, IsTF>, std::pair<FaceId, IsTF>>;

using namespace mfem;

class TotalFieldScatteredFieldSubMesher
{
public:

	TotalFieldScatteredFieldSubMesher(){};
	TotalFieldScatteredFieldSubMesher(const Mesh&, const Array<int>& marker);

	SubMesh* getTFSubMesh() { return tf_mesh_.get(); }
	SubMesh* getSFSubMesh() { return sf_mesh_.get(); }
	SubMesh* getGlobalTFSFSubMesh() { return global_submesh_.get(); }
	const SubMesh* getTFConstSubMesh() const { return tf_mesh_.get(); }
	const SubMesh* getSFConstSubMesh() const { return sf_mesh_.get(); }
	const SubMesh* getGlobalTFSFConstSubMesh() const { return global_submesh_.get(); }

private:

	void prepareSubMeshInfo(Mesh& m,   const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void setGlobalTFSFAttributesForSubMeshing(Mesh&, const Array<int>& marker);
	void storeElementToFaceInformation(const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	
	SetPairs twoPointAssignator(Mesh&, int be, bool flag);
	void assignIndividualTFSFAttsOnePoint1D(Mesh&, const Array<int>& marker);
	void assignIndividualTFSFAttsTwoPoints1D(Mesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing1D(Mesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing2D(Mesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing3D(Mesh&, const Array<int>& marker);


	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;
	std::vector<ElementId> elems_for_global_submesh_;

	std::unique_ptr<SubMesh> tf_mesh_;
	std::unique_ptr<SubMesh> sf_mesh_;
	std::unique_ptr<SubMesh> global_submesh_;

};

class NearToFarFieldSubMesher
{
public:
	NearToFarFieldSubMesher(){};
	NearToFarFieldSubMesher(const Mesh&, const FiniteElementSpace&, const Array<int>& marker);

	const SubMesh* getConstSubMesh() { return ntff_mesh_.get(); }
	SubMesh* getSubMesh() { return ntff_mesh_.get(); }
	const std::vector<El2Face> getEl2Face() { return elem_to_face_ntff_; }

private:

	void setIndividualNTFFAttributesForSubMeshing2D(Mesh& m, const Array<int>& marker);
	void setIndividualNTFFAttributesForSubMeshing3D(Mesh& m, const Array<int>& marker);
	void prepareSubMeshInfo(Mesh&, const FaceElementTransformations*, int faceId, bool el1_is_tf);
	void setAttributeForTagging(Mesh&, const FaceElementTransformations*, bool el1_is_tf);
	void storeElementToFaceInformation(const FaceElementTransformations*, int faceId, bool el1_is_tf);

	std::unique_ptr<Mesh> original_;
	std::vector<El2Face> elem_to_face_ntff_;
	std::vector<ElementId> elems_for_global_submesh_;
	std::unique_ptr<SubMesh> ntff_mesh_;

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

FaceElementTransformations* getFaceElementTransformation(Mesh& m, int be);
Vector getBarycenterOfElement(Mesh& m, int e);
Vector getBarycenterOfFaceElement(Mesh& m, int f);
Vector getNormal(FaceElementTransformations& fet);
std::pair<double, double> calculateBaryNormalProduct(Mesh& m, FaceElementTransformations& fet, int be);
double calculateCrossBaryVertexSign(Mesh& m, FaceElementTransformations& fet, int be);
double calculateFaceOrientation(Mesh& mesh, int be);
const Vector calculateTangent2D(Mesh& m, int be);
const Vector calculateNormal3D(Mesh& m, int be);
const std::pair<Vector, Vector> calculateBarycenters(Mesh& m, int be);
const Vector calculateBarycenterVector(Mesh& m, int be);

};