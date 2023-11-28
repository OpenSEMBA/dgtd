#include "Exporter.h"

namespace maxwell {

using namespace mfem;

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
	gFields_{ globalFields(global) },
	tMaps_{ TransferMaps(gFields_, fields_) }
{
	this->SetMesh(ntff_smsh_.getSubMesh());
	updateFields();
}

void NearToFarFieldDataCollection::updateFields()
{
	tMaps_.transferFields(gFields_, fields_);
}

void TransferMaps::transferFields(const globalFields& src, Fields& dst)
{
	tMapEx.Transfer(src.Ex, dst.get(E, X));
	tMapEy.Transfer(src.Ey, dst.get(E, Y));
	tMapEz.Transfer(src.Ez, dst.get(E, Z));
	tMapHx.Transfer(src.Hx, dst.get(H, X));
	tMapHy.Transfer(src.Hy, dst.get(H, Y));
	tMapHz.Transfer(src.Hz, dst.get(H, Z));
}

}