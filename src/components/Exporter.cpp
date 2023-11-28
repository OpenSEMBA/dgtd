#include "Exporter.h"

namespace maxwell {

using namespace mfem;

GlobalFields& GlobalFields::operator=(GlobalFields&& gfs)
{
	if (this == &gfs)
	{
		Ex = std::move(gfs.Ex);
		Ey = std::move(gfs.Ey);
		Ez = std::move(gfs.Ez);
		Hx = std::move(gfs.Hx);
		Hy = std::move(gfs.Hy);
		Hz = std::move(gfs.Hz);
	}
	return *this;
}

Array<int> buildSurfaceMarker(const std::vector<int>& tags, const FiniteElementSpace& fes)
{
	Array<int> res(fes.GetMesh()->bdr_attributes.Max());
	res = 0;
	for (const auto& t : tags) {
		res[t - 1] = 1;
	}
	return res;
}

void NearToFarFieldDataCollection::assignGlobalFieldsReferences(Fields& global)
{
	gFields_.Ex = global.get(E, X);
	gFields_.Ey = global.get(E, Y);
	gFields_.Ez = global.get(E, Z);
	gFields_.Hx = global.get(H, X);
	gFields_.Hy = global.get(H, Y);
	gFields_.Hz = global.get(H, Z);
}

NearToFarFieldDataCollection::NearToFarFieldDataCollection(
	const NearToFarFieldProbe& p, const DG_FECollection& fec, FiniteElementSpace& fes, Fields& global)
	: DataCollection(p.name, fes.GetMesh()),
	ntff_smsh_{ NearToFarFieldSubMesher(*fes.GetMesh(), fes, buildSurfaceMarker(p.tags, fes)) },
	sfes_{ std::make_unique<FiniteElementSpace>(ntff_smsh_.getSubMesh(), &fec) },
	fields_{ Fields(*sfes_.get()) },
	gFields_{ GlobalFields(global) },
	tMaps_{ TransferMaps(gFields_, fields_) }
{
	this->SetMesh(ntff_smsh_.getSubMesh());
	updateFields();
}

NearToFarFieldDataCollection& NearToFarFieldDataCollection::operator=(NearToFarFieldDataCollection&& dc)
{
	if (this != &dc)
	{
		ntff_smsh_ = std::move(dc.ntff_smsh_);
		sfes_ = std::move(dc.sfes_);
		fields_ = std::move(dc.fields_);
		gFields_ = std::move(dc.gFields_);
		tMaps_ = std::move(dc.tMaps_);
	}

	return *this;

}

void NearToFarFieldDataCollection::updateFields()
{
	tMaps_.transferFields(gFields_, fields_);
}

void TransferMaps::transferFields(const GlobalFields& src, Fields& dst)
{
	tMapEx.Transfer(src.Ex, dst.get(E, X));
	tMapEy.Transfer(src.Ey, dst.get(E, Y));
	tMapEz.Transfer(src.Ez, dst.get(E, Z));
	tMapHx.Transfer(src.Hx, dst.get(H, X));
	tMapHy.Transfer(src.Hy, dst.get(H, Y));
	tMapHz.Transfer(src.Hz, dst.get(H, Z));
}

}