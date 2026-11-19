#include "Model.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include "PMLProfiles.h"

#include <memory>
#include <string>

namespace maxwell {

using namespace mfem;


std::map<GlobalElementId, Position> buildSerialElem2CenterMap(Mesh& mesh){
	std::map<GlobalElementId, Position> res;
	for (auto e = 0; e < mesh.GetNE(); e++){
		Vector center;
		mesh.GetElementCenter(e, center);
		res[e] = center;
	}
	return res;
}


std::map<LocalElementId, Position> buildPartitionElem2CenterMap(ParMesh& pmesh){
	std::map<LocalElementId, Position> res;
	for (auto e = 0; e < pmesh.GetNE(); e++){
		Vector center;
		pmesh.GetElementCenter(e, center);
		res[e] = center;
	}
	return res;
}

std::map<GlobalElementId, LocalElementId> buildGlobalToPartitionLocalElementMap(const std::map<GlobalElementId, Position>& serial, const std::map<LocalElementId, Position>& local)
{
	double tol = 1e-5;
	std::map<GlobalElementId, LocalElementId> res;
	for (const auto& [glob_el_id, center_to_find] : serial)	{
		for (const auto& [loc_el_id, local_cent] : local){
			if( center_to_find.DistanceTo(local_cent) <= tol){
				res[glob_el_id] = loc_el_id;
			}
		}
	}
	return res;
}

void ensureElementTypeIsSame(const Mesh& mesh)
{
	const auto type = mesh.GetElementType(0);
	for (auto e = 1; e < mesh.GetNE(); e++){
		if (type != mesh.GetElementType(e)){
			throw std::runtime_error("Parallelization only works if all mesh elements are of the same type");
		}	
	}
}

Model::Model(Mesh& mesh, const GeomTagToMaterialInfo& matInfo, const GeomTagToBoundaryInfo& bdrInfo, int* partitioning, MPI_Comm comm)
{

	serialMesh_ = Mesh(mesh);
	ensureElementTypeIsSame(mesh);

	pmesh_ = ParMesh(comm, serialMesh_, partitioning);

	if (matInfo.gt2m.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0, 0.0));
	}
	else {
		attToMatMap_ = matInfo.gt2m;
	}

	if (matInfo.gt2bm.size() != 0){
		attToBdrMatMap_ = matInfo.gt2bm;
	}

	Array<int> f2bdr;
	if (partitioning != nullptr){
		f2bdr = pmesh_.GetFaceToBdrElMap();
	}
	else{
		f2bdr = serialMesh_.GetFaceToBdrElMap();
	}

	for (auto i = bdrInfo.gt2b.begin(); i != bdrInfo.gt2b.end(); i++) {
		faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
	}
	attToBdrMap_ = bdrInfo.gt2b;
	
	if (bdrInfo.gt2ib.size() != 0)
	{
		for (auto i = bdrInfo.gt2ib.begin(); i != bdrInfo.gt2ib.end(); i++){
			if (i->second == BdrCond::PEC || i->second == BdrCond::PMC || i->second == BdrCond::SMA || i->second == BdrCond::SGBC) {
				faceToGeomTag_.insert(std::make_pair(f2bdr.Find(i->first - 1), i->first));
				attToIntBdrMap_.insert(std::make_pair( i->first, i->second ));
			}
		}
	}
	attToIntBdrMap_ = bdrInfo.gt2ib;

	assembleGeomTagToTypeMap(attToBdrMap_, false);
	assembleGeomTagToTypeMap(attToIntBdrMap_, true);
	assembleBdrToMarkerMaps();

	if (partitioning != nullptr) {
		const int ne = serialMesh_.GetNE();
		element_partition_.SetSize(ne);
		for (int i = 0; i < ne; ++i) {
			element_partition_[i] = partitioning[i];
		}
	}

}

void Model::setPMLProperties(const std::vector<PMLProperties>& in)
{
	pml_props_ = in;
	rejectDispersiveOnPML();
}

void Model::setDebyeProperties(const std::vector<DebyeProperties>& in)
{
	debye_ = in;
	debye_by_tag_.clear();
	for (const DebyeProperties& pole : debye_) {
		debye_by_tag_.insert_or_assign(pole.geom_tag, pole);
	}
	rejectDispersiveOnPML();
}

const DebyeProperties* Model::findDebye(Attribute tag) const
{
	const auto it = debye_by_tag_.find(tag);
	if (it == debye_by_tag_.end()) {
		return nullptr;
	}
	return &it->second;
}

void Model::setLorentzProperties(const std::vector<LorentzProperties>& in)
{
	lorentz_ = in;
	if (serialMesh_.GetNE() > 0) {
		for (LorentzProperties& pole : lorentz_) {
			if (pole.hasRegionBounds()) {
				continue;
			}
			double min_x = std::numeric_limits<double>::infinity();
			double max_x = -std::numeric_limits<double>::infinity();
			mfem::Array<int> verts;
			for (int el = 0; el < serialMesh_.GetNE(); ++el) {
				if (serialMesh_.GetAttribute(el) != pole.geom_tag) {
					continue;
				}
				serialMesh_.GetElementVertices(el, verts);
				for (int v = 0; v < verts.Size(); ++v) {
					const double vx = serialMesh_.GetVertex(verts[v])[0];
					min_x = std::min(min_x, vx);
					max_x = std::max(max_x, vx);
				}
			}
			if (max_x > min_x) {
				pole.region_min_x = min_x;
				pole.region_max_x = max_x;
			}
		}
	}
	lorentz_by_tag_.clear();
	for (const LorentzProperties& pole : lorentz_) {
		lorentz_by_tag_.insert_or_assign(pole.geom_tag, pole);
	}
	rejectDispersiveOnPML();
}

namespace {

long long countGridFunctionValues(const std::string& file)
{
	std::ifstream input(file);
	if (!input) {
		return -1;
	}
	std::string line;
	bool in_values = false;
	long long count = 0;
	while (std::getline(input, line)) {
		if (!in_values) {
			if (line.rfind("Ordering", 0) == 0) {
				in_values = true;
			}
			continue;
		}
		std::istringstream row(line);
		std::string token;
		while (row >> token) {
			++count;
		}
	}
	return count;
}

std::shared_ptr<mfem::ParGridFunction> loadDistributedGridFunction(
	const std::string& file, mfem::Mesh& serial_mesh, mfem::ParMesh& pmesh,
	mfem::ParFiniteElementSpace& pfes)
{
	const long long values = countGridFunctionValues(file);
	if (values < 0) {
		throw std::runtime_error(
			"Lorentz grid function file cannot be opened: " + file);
	}

	std::ifstream input(file);
	mfem::GridFunction serial(&serial_mesh, input);
	const long long expected = serial.FESpace()->GetVSize();
	if (values != expected) {
		throw std::runtime_error(
			"Lorentz grid function " + file + " has " + std::to_string(values) +
			" values, but the case mesh and FES need " + std::to_string(expected) +
			". Regenerate it with plasma_kit for this mesh, order, basis type and refinement.");
	}
	if (serial.FESpace()->GetVDim() != 1) {
		throw std::runtime_error(
			"Lorentz grid function file must be scalar: " + file);
	}
	if (serial.FESpace()->GetNE() != serial_mesh.GetNE()) {
		throw std::runtime_error(
			"Lorentz grid function " + file + " was built on a mesh with " +
			std::to_string(serial.FESpace()->GetNE()) + " elements, but the case mesh has " +
			std::to_string(serial_mesh.GetNE()) +
			". Regenerate it with plasma_kit for this mesh and refinement.");
	}
	if (serial.FESpace()->GetFE(0)->GetOrder() != pfes.GetFE(0)->GetOrder()) {
		throw std::runtime_error(
			"Lorentz grid function " + file + " has order " +
			std::to_string(serial.FESpace()->GetFE(0)->GetOrder()) + " but the case uses order " +
			std::to_string(pfes.GetFE(0)->GetOrder()) + ".");
	}
	if (serial.FESpace()->GetFE(0)->GetDof() != pfes.GetFE(0)->GetDof()) {
		throw std::runtime_error(
			"Lorentz grid function " + file +
			" has a different number of DoFs per element than the case FES.");
	}

	auto parallel = std::make_shared<mfem::ParGridFunction>(&pfes);
	mfem::Array<int> serial_dofs, local_dofs;
	for (int e = 0; e < pmesh.GetNE(); ++e) {
		const int global = static_cast<int>(pmesh.GetGlobalElementNum(e));
		if (global < 0 || global >= serial_mesh.GetNE()) {
			throw std::runtime_error(
				"Lorentz grid function " + file + ": invalid global element index.");
		}
		serial.FESpace()->GetElementDofs(global, serial_dofs);
		pfes.GetElementDofs(e, local_dofs);
		if (serial_dofs.Size() != local_dofs.Size()) {
			throw std::runtime_error(
				"Lorentz grid function " + file +
				": element DoF count mismatch during distribution.");
		}
		for (int i = 0; i < local_dofs.Size(); ++i) {
			(*parallel)(local_dofs[i]) = serial(serial_dofs[i]);
		}
	}
	return parallel;
}

double maxValue(const mfem::ParGridFunction& gf)
{
	double value = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < gf.Size(); ++i) {
		value = std::max(value, gf(i));
	}
	return value;
}

double minValue(const mfem::ParGridFunction& gf)
{
	double value = std::numeric_limits<double>::infinity();
	for (int i = 0; i < gf.Size(); ++i) {
		value = std::min(value, gf(i));
	}
	return value;
}

} // namespace

void Model::loadLorentzGridFunctions(mfem::ParFiniteElementSpace& fes)
{
	for (auto& [tag, pole] : lorentz_by_tag_) {
		(void)tag;
		if (!pole.omega_p_file.empty()) {
			pole.omega_p_gf = loadDistributedGridFunction(
				pole.omega_p_file, serialMesh_, pmesh_, fes);
			if (minValue(*pole.omega_p_gf) < 0.0) {
				throw std::runtime_error(
					"Lorentz omega_p grid function has negative values: " +
					pole.omega_p_file);
			}
			pole.omega_p = maxValue(*pole.omega_p_gf);
		}
		if (!pole.omega_1_file.empty()) {
			pole.omega_1_gf = loadDistributedGridFunction(
				pole.omega_1_file, serialMesh_, pmesh_, fes);
			pole.omega_1 = maxValue(*pole.omega_1_gf);
		}
		if (!pole.gamma_file.empty()) {
			pole.gamma_gf = loadDistributedGridFunction(
				pole.gamma_file, serialMesh_, pmesh_, fes);
			if (minValue(*pole.gamma_gf) < 0.0) {
				throw std::runtime_error(
					"Lorentz gamma grid function has negative values: " +
					pole.gamma_file);
			}
			pole.gamma = maxValue(*pole.gamma_gf);
		}
	}
	for (LorentzProperties& pole : lorentz_) {
		const auto it = lorentz_by_tag_.find(pole.geom_tag);
		if (it != lorentz_by_tag_.end()) {
			pole.omega_p_gf = it->second.omega_p_gf;
			pole.omega_1_gf = it->second.omega_1_gf;
			pole.gamma_gf = it->second.gamma_gf;
			pole.omega_p = it->second.omega_p;
			pole.omega_1 = it->second.omega_1;
			pole.gamma = it->second.gamma;
		}
	}
}
void Model::rejectDispersiveOnPML() const
{
	if (pml_props_.empty() || (debye_.empty() && lorentz_.empty())) {
		return;
	}
	for (const PMLProperties& props : pml_props_) {
		for (const GeomTag tag : props.geom_tags) {
			if (debye_by_tag_.count(tag) != 0 || lorentz_by_tag_.count(tag) != 0) {
				throw std::runtime_error(kDispersiveOnPmlNotAllowed);
			}
		}
	}
}

const LorentzProperties* Model::findLorentz(Attribute tag) const
{
	const auto it = lorentz_by_tag_.find(tag);
	if (it == lorentz_by_tag_.end()) {
		return nullptr;
	}
	return &it->second;
}

void Model::assembleBdrToMarkerMaps()
{
	const std::pair<BdrCond, const mfem::Array<int>&> bdrEntries[] = {
		{BdrCond::PEC, pecMarker_}, {BdrCond::PMC, pmcMarker_},
		{BdrCond::SMA, smaMarker_},
		{BdrCond::SGBC, sgbc_Marker_}
	};
	for (const auto& [cond, marker] : bdrEntries) {
		if (marker.Size() != 0) bdrToMarkerMap_.insert({cond, marker});
	}

	const std::pair<BdrCond, const mfem::Array<int>&> intBdrEntries[] = {
		{BdrCond::PEC, intpecMarker_}, {BdrCond::PMC, intpmcMarker_},
		{BdrCond::SMA, intsmaMarker_},
		{BdrCond::SGBC, intsgbc_Marker_}
	};
	for (const auto& [cond, marker] : intBdrEntries) {
		if (marker.Size() != 0) intBdrToMarkerMap_.insert({cond, marker});
	}
}

std::size_t Model::numberOfMaterials() const
{
	return attToMatMap_.size();
}

std::size_t Model::numberOfBoundaryMaterials() const
{
	return attToBdrMap_.size();
}

mfem::Vector Model::initialiseGeomTagVector() const
{
	int size{ 0 };
	for (auto const& [geomTag, mat] : attToMatMap_) {
		if (geomTag >= size) {
			size = geomTag;
		}
	}
	assert(size > 0);
	mfem::Vector res(size);
	res = 0.0;
	return res;
}

mfem::Vector Model::buildEpsMuPiecewiseVector(const FieldType& f) const
{
	auto res{ initialiseGeomTagVector() };

	for (auto const& [geomTag, mat] : attToMatMap_) {
		if (geomTag <= 0 || geomTag > res.Size()) {
			throw std::runtime_error(
				"Material tag " + std::to_string(geomTag)
				+ " is outside 1.." + std::to_string(res.Size()) + ".");
		}
		switch (f) {
		case FieldType::E:
			res[geomTag - 1] = mat.getPermittivity();
			break;
		case FieldType::H:
			res[geomTag - 1] = mat.getPermeability();
			break;
		}
	}

	return res;
}

mfem::Vector Model::buildSigmaPiecewiseVector() const
{
	auto res{ initialiseGeomTagVector() };

	for (auto const& [geomTag, mat] : attToMatMap_) {
		if (geomTag <= 0 || geomTag > res.Size()) {
			throw std::runtime_error(
				"Material tag " + std::to_string(geomTag)
				+ " is outside 1.." + std::to_string(res.Size()) + ".");
		}
		res[geomTag - 1] = mat.getConductivity();
	}

	return res;
}

void initMarker(BoundaryMarker& marker, const int size)
{
	marker.SetSize(size);
	marker = 0;
}

void Model::assembleGeomTagToTypeMap(
	std::map<GeomTag, BdrCond>& geomTagToCond,
	bool isInterior)
{
	for (const auto& [geomTag, bdr] : geomTagToCond) {

		if (geomTag <= 0) {
			throw std::runtime_error("geomTag <= 0 in GeomTagToTypeMap assembly.");
		}

		auto& marker{ getMarker(bdr, isInterior) };

		if (marker.Size() == 0) {
			const int n_attr = pmesh_.bdr_attributes.Size() == 0
				? 0 : pmesh_.bdr_attributes.Max();
			if (n_attr <= 0) {
				throw std::runtime_error(
					"Boundary tag " + std::to_string(geomTag)
					+ " is set but the mesh has no boundary attributes.");
			}
			initMarker(getMarker(bdr, isInterior), n_attr);
		}

		if (geomTag > marker.Size()) {
			throw std::runtime_error(
				"Boundary tag " + std::to_string(geomTag)
				+ " is outside the mesh boundary-attribute range 1.."
				+ std::to_string(marker.Size()) + ".");
		}

		marker[geomTag - 1] = 1;
	}
}

BoundaryMarker& Model::getMarker(const BdrCond& bdrCond, bool isInterior)
{
	switch (bdrCond) {
	case BdrCond::PEC:
		switch (isInterior) {
			case true:
				return intpecMarker_;
			case false:
				return pecMarker_;
		}
		break;
	case BdrCond::PMC:
		switch (isInterior) {
			case true:
				return intpmcMarker_;
			case false:
				return pmcMarker_;
		}
		break;
	case BdrCond::SMA:
		switch (isInterior) {
			case true:
				return intsmaMarker_;
			case false:
				return smaMarker_;
		}
		break;
	case BdrCond::TotalFieldIn:
		return tfsfMarker_;
		break;
	case BdrCond::SGBC:
		switch (isInterior) {
			case true:
				return intsgbc_Marker_;
			case false:
				return sgbc_Marker_;
		}
		break;
	default:
		throw std::runtime_error("Wrong BdrCond in getMarkerForBdrCond.");
	}
}

bool Model::isPMLAttribute(GeomTag tag) const
{
	return getPMLPropertiesForTag(tag) != nullptr;
}

const PMLProperties* Model::getPMLPropertiesForTag(GeomTag tag) const
{
	for (const auto& props : pml_props_) {
		for (const GeomTag t : props.geom_tags) {
			if (t == tag) {
				return &props;
			}
		}
	}
	return nullptr;
}

void Model::initializePMLProfiles(int mpi_rank, int fe_order)
{
	if (pml_props_.empty()) {
		return;
	}
	// Build on the full serial mesh so vacuum–PML interface coords and
	// region_max_depth (L) are global — not truncated to a partition. σ(x)
	// evaluation during ParBilinearForm assembly uses Attribute + physical
	// coordinates (see PMLProfileData::evaluateAtTransform), not ElementNo.
	pml_profiles_ = std::make_shared<PMLProfileData>(serialMesh_, pml_props_, fe_order);
	pml_profiles_->printDiagnostics(mpi_rank);
}

mfem::Array<int> Model::buildPMLVolumeMarker() const
{
	const int max_attr = serialMesh_.attributes.Max();
	mfem::Array<int> marker(max_attr);
	marker = 0;
	for (const auto& props : pml_props_) {
		for (const GeomTag tag : props.geom_tags) {
			if (tag <= 0 || tag > max_attr) {
				throw std::runtime_error("PML material tag is out of mesh attribute range.");
			}
			marker[tag - 1] = 1;
		}
	}
	return marker;
}

}
