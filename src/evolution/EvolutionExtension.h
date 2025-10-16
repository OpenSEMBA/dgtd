#include "components/Model.h"
#include "HesthavenEvolutionMethods.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

struct SBCProperties{

    size_t num_of_segments = 10;
    size_t order = 1;
    size_t material_width = 1e-4;

    SBCProperties(size_t segnum, size_t o, size_t mat_w) : 
    num_of_segments(segnum), order(o), material_width(mat_w){}

};

class SBCManager{
public:

SBCManager(Model&, FiniteElementSpace&, const SBCProperties&);

private:

void findDoFPairs(Model&, FiniteElementSpace&);

std::unique_ptr<Mesh> mesh_;
std::unique_ptr<DG_FECollection> fec_;
std::unique_ptr<FiniteElementSpace> fes_;

std::vector<std::pair<NodeId, NodeId>> dof_pairs_;

};

}