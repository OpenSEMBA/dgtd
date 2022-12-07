#pragma once 

#include "geometry/mesh/Unstructured.h"

namespace SEMBA {
namespace dgtd {
namespace Mesh {

class Volume : public SEMBA::Geometry::Mesh::Unstructured {
public:
    Volume(const Volume& meshVol);
    Volume(const Unstructured& uns);
    
    Volume& operator=(const Volume& param);
    
    vector<vector<Geometry::ElemId>> getPartitionsIds(
            const size_t nDivisions,
            const vector<pair<Geometry::ElemId,int>> idWeights =
                    vector<pair<Geometry::ElemId,int>>(),
            const Math::Real* taskPower = NULL) const;

    const Geometry::Graph::Connectivities* getConnectivities() const;
private:
    Geometry::Graph::Connectivities* connectivities_;
};

} 
}
}

#endif
