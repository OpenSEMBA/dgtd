#include "Exporter.h"

void maxwell::NearToFarFieldDataCollection::assignGlobalFieldsReferences(Fields& global)
{
	gFields_.Ex = global.get(E, X);
	gFields_.Ey = global.get(E, Y);
	gFields_.Ez = global.get(E, Z);
	gFields_.Hx = global.get(H, X);
	gFields_.Hy = global.get(H, Y);
	gFields_.Hz = global.get(H, Z);
}

maxwell::NearToFarFieldDataCollection::NearToFarFieldDataCollection(const std::string& name, mfem::FiniteElementSpace& subFes, Fields& global) 
	: DataCollection(name, subFes.GetMesh()),
	sfes_(subFes),
	fields_{Fields(sfes_)},
	tMaps_{ TransferMaps(global, fields_) },
	gFields_ { globalFields(global) }
{
	updateFields();
}

void maxwell::NearToFarFieldDataCollection::updateFields()
{
	tMaps_.transferFields(gFields_, fields_);
}

void maxwell::TransferMaps::transferFields(const globalFields& src, Fields& dst)
{
	Ex.Transfer(src.Ex, dst.get(E, X));
	Ey.Transfer(src.Ey, dst.get(E, Y));
	Ez.Transfer(src.Ez, dst.get(E, Z));
	Hx.Transfer(src.Hx, dst.get(H, X));
	Hy.Transfer(src.Hy, dst.get(H, Y));
	Hz.Transfer(src.Hz, dst.get(H, Z));
}
