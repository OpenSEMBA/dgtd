#pragma once

#include <evolution/Fields.h>
#include <components/Types.h>

namespace maxwell {

struct globalFields {

	mfem::GridFunction& Ex;
	mfem::GridFunction& Ey;
	mfem::GridFunction& Ez;
	mfem::GridFunction& Hx;
	mfem::GridFunction& Hy;
	mfem::GridFunction& Hz;

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

	mfem::TransferMap Ex;
	mfem::TransferMap Ey;
	mfem::TransferMap Ez;
	mfem::TransferMap Hx;
	mfem::TransferMap Hy;
	mfem::TransferMap Hz;

	TransferMaps(Fields& src, Fields& dst) :
		Ex{ mfem::TransferMap(src.get(E, X), dst.get(E, X)) },
		Ey{ mfem::TransferMap(src.get(E, Y), dst.get(E, Y)) },
		Ez{ mfem::TransferMap(src.get(E, Z), dst.get(E, Z)) },
		Hx{ mfem::TransferMap(src.get(H, X), dst.get(H, X)) },
		Hy{ mfem::TransferMap(src.get(H, Y), dst.get(H, Y)) },
		Hz{ mfem::TransferMap(src.get(H, Z), dst.get(H, Z)) } 
	{}

	void transferFields(const globalFields&, Fields&);
};

class NearToFarFieldDataCollection : public mfem::DataCollection
{
public:

	NearToFarFieldDataCollection(const std::string&, mfem::FiniteElementSpace&, Fields&);

	mfem::GridFunction& getCollectionField(const FieldType& f, const Direction& d)  { return fields_.get(f, d) ; }
	void updateFields();

private:

	void assignGlobalFieldsReferences(Fields& global);
	
	mfem::FiniteElementSpace sfes_;
	Fields fields_;
	TransferMaps tMaps_;
	globalFields gFields_;

};

}