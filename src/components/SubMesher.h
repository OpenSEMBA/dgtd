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
	TotalFieldScatteredFieldSubMesher(const ParMesh&, const Array<int>& marker);

	ParSubMesh* getTFSubMesh() { return tf_mesh_.get(); }
	ParSubMesh* getSFSubMesh() { return sf_mesh_.get(); }
	ParSubMesh* getGlobalTFSFSubMesh() { return global_submesh_.get(); }
	const ParSubMesh* getTFConstSubMesh() const { return tf_mesh_.get(); }
	const ParSubMesh* getSFConstSubMesh() const { return sf_mesh_.get(); }
	const ParSubMesh* getGlobalTFSFConstSubMesh() const { return global_submesh_.get(); }

private:

	void prepareSubMeshInfo(ParMesh& m,   const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	void setAttributeForTagging(ParMesh&, const FaceElementTransformations*, bool el1_is_tf);
	void setGlobalTFSFAttributesForSubMeshing(ParMesh&, const Array<int>& marker);
	void storeElementToFaceInformation(const FaceElementTransformations*, const std::pair<int, int> facesId, bool el1_is_tf);
	
	SetPairs twoPointAssignator(ParMesh&, int be, bool flag);
	void assignIndividualTFSFAttsOnePoint1D(ParMesh&, const Array<int>& marker);
	void assignIndividualTFSFAttsTwoPoints1D(ParMesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing1D(ParMesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing2D(ParMesh&, const Array<int>& marker);
	void setIndividualTFSFAttributesForSubMeshing3D(ParMesh&, const Array<int>& marker);


	std::vector<El2Face> elem_to_face_tf_;
	std::vector<El2Face> elem_to_face_sf_;
	std::vector<ElementId> elems_for_global_submesh_;

	std::unique_ptr<ParSubMesh> tf_mesh_;
	std::unique_ptr<ParSubMesh> sf_mesh_;
	std::unique_ptr<ParSubMesh> global_submesh_;

};

class NearToFarFieldSubMesher
{
public:
	NearToFarFieldSubMesher(){};
	NearToFarFieldSubMesher(const Mesh&, const ParFiniteElementSpace&, const Array<int>& marker);

	const ParSubMesh* getConstSubMesh() { return ntff_mesh_.get(); }
	ParSubMesh* getSubMesh() { return ntff_mesh_.get(); }
	const std::vector<El2Face> getEl2Face() { return elem_to_face_ntff_; }

private:

	void setSurfaceAttributesForSubMesh2D(ParMesh& m, const Array<int>& marker);
	void setSurfaceAttributesForSubMesh3D(ParMesh& m, const Array<int>& marker);
	void prepareSubMeshInfo(ParMesh&, const FaceElementTransformations*, int faceId, bool el1_is_tf);
	void setAttributeForTagging(ParMesh&, const FaceElementTransformations*, bool el1_is_tf);
	void storeElementToFaceInformation(const FaceElementTransformations*, int faceId, bool el1_is_tf);

	std::unique_ptr<ParMesh> original_;
	std::vector<El2Face> elem_to_face_ntff_;
	std::vector<ElementId> elems_for_global_submesh_;
	std::unique_ptr<ParSubMesh> ntff_mesh_;

};

class MaxwellTransferMap
{
public:
	
	MaxwellTransferMap(const ParGridFunction& src, const ParGridFunction& dst);
	
	void TransferAdd(const ParGridFunction& src, ParGridFunction& dst) const;
	void TransferSub(const ParGridFunction& src, ParGridFunction& dst) const;

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
double buildFaceOrientation(Mesh& mesh, int be);
const Vector buildTangent2D(Mesh& m, int be);
const Vector buildNormal3D(Mesh& m, int be);
const std::pair<Vector, Vector> calculateBarycenters(Mesh& m, int be);
const Vector buildElem1ToElem2BarycenterVector(Mesh& m, int be);

};