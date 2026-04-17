#include "components/RCSManager.h"
#include "components/FarField.h"
#include <string>
#include "components/Types.h"
#include "FieldUtils.h"


namespace maxwell{

using namespace mfem;

using Frequency = double;

struct ComplexField {
    std::unique_ptr<mfem::GridFunction> real;
    std::unique_ptr<mfem::GridFunction> imag;
};

struct FrequencyFields {
    ComplexField Ex, Ey, Ez, Hx, Hy, Hz;
};

using FrequencyToFields = std::map<Frequency, FrequencyFields>;

struct CaseInfo {
    const std::string data_path;
    const std::string json_path;
    std::vector<std::string> rank_paths;
    
    CaseInfo(const std::string& d, const std::string& j);
};

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