#pragma once

#include "Ordering.h"
#include "mesh/Volume.h"

namespace SEMBA {
namespace dgtd {
namespace Cell {

class Group {
public:
	vector<CellTet<ORDER_N>*> cell;
	vector<CellTet4<ORDER_N> > linTet;
	vector<CellTet10<ORDER_N> > quadTet;

	Group(const VolumeMesh& mesh, const PMGroup& pMGroup);
	
	const CellTet<ORDER_N>* operator()(const size_t i) const;
	const CellTet<ORDER_N>* getPtrToCell(const Geometry::VolR* elem) const;
	const CellTet<ORDER_N>* getPtrToCellWithId(const Geometry::ElemId&) const;

private:
	size_t cellOffsetId;
	void buildNodalMaps(const Geometry::Graph::Connectivities& map);
	void check(const Geometry::Graph::Connectivities& map) const;
	void checkNodalMaps(const Geometry::Graph::Connectivities& map) const;
};

}
}
}