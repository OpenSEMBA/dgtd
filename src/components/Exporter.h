#pragma once

#include <evolution/Fields.h>
#include <components/Types.h>
#include <components/SubMesher.h>
#include <components/Probes.h>

namespace maxwell {

using namespace mfem;

struct GlobalFields {

	GridFunction& Ex;
	GridFunction& Ey;
	GridFunction& Ez;
	GridFunction& Hx;
	GridFunction& Hy;
	GridFunction& Hz;

	GlobalFields(Fields& global) :
		Ex{ global.get(E, X) },
		Ey{ global.get(E, Y) },
		Ez{ global.get(E, Z) },
		Hx{ global.get(H, X) },
		Hy{ global.get(H, Y) },
		Hz{ global.get(H, Z) }
	{}

	GlobalFields(const GlobalFields&) = delete;
	GlobalFields(GlobalFields&&) = default;

	GlobalFields& operator=(const GlobalFields&) = delete;
	GlobalFields& operator=(GlobalFields&& gfs);

	~GlobalFields() = default;


};

struct TransferMaps {

	TransferMap tMapEx;
	TransferMap tMapEy;
	TransferMap tMapEz;
	TransferMap tMapHx;
	TransferMap tMapHy;
	TransferMap tMapHz;

	TransferMaps(GlobalFields& src, Fields& dst) :
		tMapEx{ TransferMap(src.Ex, dst.get(E, X)) },
		tMapEy{ TransferMap(src.Ey, dst.get(E, Y)) },
		tMapEz{ TransferMap(src.Ez, dst.get(E, Z)) },
		tMapHx{ TransferMap(src.Hx, dst.get(H, X)) },
		tMapHy{ TransferMap(src.Hy, dst.get(H, Y)) },
		tMapHz{ TransferMap(src.Hz, dst.get(H, Z)) }
	{}

	void transferFields(const GlobalFields&, Fields&);
};

class NearToFarFieldDataCollection : public DataCollection
{
public:

	NearToFarFieldDataCollection(const NearToFarFieldProbe&, const DG_FECollection& fec, FiniteElementSpace& fes, Fields&);

	NearToFarFieldDataCollection(const NearToFarFieldDataCollection&) = delete;
	NearToFarFieldDataCollection(NearToFarFieldDataCollection&&);

	NearToFarFieldDataCollection& operator=(const NearToFarFieldDataCollection&) = delete;
	NearToFarFieldDataCollection& operator=(NearToFarFieldDataCollection&&);

	~NearToFarFieldDataCollection() = default;

	GridFunction& getCollectionField(const FieldType& f, const Direction& d)  { return fields_.get(f, d) ; }
	void updateFields();

private:

	void assignGlobalFieldsReferences(Fields& global);
	
	NearToFarFieldSubMesher ntff_smsh_;
	std::unique_ptr<FiniteElementSpace> sfes_;
	Fields fields_;
	GlobalFields gFields_;
	TransferMaps tMaps_;

};

}