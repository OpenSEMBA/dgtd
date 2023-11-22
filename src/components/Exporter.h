#pragma once

#include <evolution/Fields.h>
#include <components/Types.h>
#include <components/SubMesher.h>
#include <components/Probes.h>

namespace maxwell {

using namespace mfem;

struct globalFields {

	GridFunction& Ex;
	GridFunction& Ey;
	GridFunction& Ez;
	GridFunction& Hx;
	GridFunction& Hy;
	GridFunction& Hz;

	globalFields(Fields& global) :
		Ex{ global.get(E,X) },
		Ey{ global.get(E,X) },
		Ez{ global.get(E,X) },
		Hx{ global.get(E,X) },
		Hy{ global.get(E,X) },
		Hz{ global.get(E,X) }
	{}
};

struct TransferMaps {

	TransferMap tMapEx;
	TransferMap tMapEy;
	TransferMap tMapEz;
	TransferMap tMapHx;
	TransferMap tMapHy;
	TransferMap tMapHz;

	TransferMaps(globalFields& src, Fields& dst) :
		tMapEx{ TransferMap(src.Ex, dst.get(E, X)) },
		tMapEy{ TransferMap(src.Ey, dst.get(E, Y)) },
		tMapEz{ TransferMap(src.Ez, dst.get(E, Z)) },
		tMapHx{ TransferMap(src.Hx, dst.get(H, X)) },
		tMapHy{ TransferMap(src.Hy, dst.get(H, Y)) },
		tMapHz{ TransferMap(src.Hz, dst.get(H, Z)) }
	{}

	void transferFields(const globalFields&, Fields&);
};

class NearToFarFieldDataCollection : public DataCollection
{
public:

	NearToFarFieldDataCollection(const NearToFarFieldProbe&, DG_FECollection& fec, FiniteElementSpace& fes, Fields&);

	GridFunction& getCollectionField(const FieldType& f, const Direction& d)  { return fields_.get(f, d) ; }
	void updateFields();

private:

	void assignGlobalFieldsReferences(Fields& global);
	
	NearToFarFieldSubMesher ntff_smsh_;
	FiniteElementSpace sfes_;
	Fields fields_;
	globalFields gFields_;
	TransferMaps tMaps_;

};

}