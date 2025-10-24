#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

GeomTagToMaterial getSBCSolverGeomTagToMaterialFromGlobal(Model& global_model);

using NbrPairs = std::pair<double, double>;

struct FluxRows
{
    Eigen::VectorXcd row_left_first;
    Eigen::VectorXcd row_left_second;
    Eigen::VectorXcd row_right_first;
    Eigen::VectorXcd row_right_second;
};

using FieldComponentToFluxRows = std::map<std::pair<FieldType,Direction>,FluxRows>;
using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

class SBCSolver{
public:

    SBCSolver(Model&, ParFiniteElementSpace&, const SBCProperties&);

    void setTargetTime(Time& t) { target_time = t; }
    void setPreTime(Time& t) { pre_time = t; }
    void assignGlobalFields(const Fields<ParFiniteElementSpace,ParGridFunction>* g_fields);
    
private:
    
    SBCProperties sbcp_;
    
    std::vector<std::pair<NodeId, NodeId>> dof_pairs_;
    
    Model model_;
    
    Time target_time;
    Time pre_time = 0.0;
    Time dt_;
    
    SolverOptions opts_;

    FieldComponentToFluxRows nodal_to_modal_rows_;
    FieldComponentToFluxRows modal_to_nodal_rows_;

    ModalValues modal_values_;
    NodalValues nodal_values_;

    std::vector<NodeId> target_ids_;
    
    const Fields<ParFiniteElementSpace, ParGridFunction>* global_fields_;
    
    void findDoFPairs(Model&, ParFiniteElementSpace&);

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);
Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnGlobalOperator(const SparseMatrix& mat);

}