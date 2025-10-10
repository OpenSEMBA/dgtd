#include "components/RCSManager.h"
#include "components/FarField.h"
#include <string>
#include "components/Types.h"


namespace maxwell{

using namespace mfem;

struct ComplexField {
    std::unique_ptr<mfem::GridFunction> re;
    std::unique_ptr<mfem::GridFunction> im;
};

struct FrequencyFields {
    ComplexField Ex, Ey, Ez, Hx, Hy, Hz;
};

struct CaseInfo {
    const std::string data_path;
    const std::string json_path;
    std::vector<std::string> rank_paths;
    
    CaseInfo(const std::string& d, const std::string& j);
};

struct TimeFields {
    std::unique_ptr<GridFunction> Ex;
    std::unique_ptr<GridFunction> Ey;
    std::unique_ptr<GridFunction> Ez;
    std::unique_ptr<GridFunction> Hx;
    std::unique_ptr<GridFunction> Hy;
    std::unique_ptr<GridFunction> Hz;
};

using TimeToFields = std::map<double, TimeFields>;

class FieldSuperposition {
public:

    FieldSuperposition(const CaseInfo& c1, const CaseInfo& c2, const CaseInfo& c_compare, const Frequency);

private:

    void loadRankPaths();
    Mesh loadCaseMesh(const size_t& rank);

    CaseInfo c1_, c2_, c_compare_;
    TimeToFields t2f1_, t2f2_, t2f_compare_;

};

}