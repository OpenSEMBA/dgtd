#include "components/RCSManager.h"
#include "components/FarField.h"
#include <string>
#include "components/Types.h"


namespace maxwell{


struct CaseInfo {
    const std::string data_path;
    const std::string json_path;
    std::vector<std::string> rank_paths;

    CaseInfo(const std::string& d, const std::string& j);
};

struct TimeFields {
    GridFunction Ex, Ey, Ez,
                 Hx, Hy, Hz;
};

class FieldSuperposition {
public:

    FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare);

private:

    void loadRankPaths();
    Mesh loadCaseMesh(const size_t& rank);

    CaseInfo c1_, c2_, c_compare_;
    std::map<Time, TimeFields> t2f1_, t2f2_, t2f_compare_;

};

}