#pragma once

#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

GeomTagToMaterial getSBCSolverGeomTagToMaterialFromGlobal(Model& global_model);

using GlobalId = NodeId;
using LocalId = NodeId;

struct FluxRows
{
    Eigen::VectorXcd row_left_first;
    Eigen::VectorXcd row_left_second;
    Eigen::VectorXcd row_right_first;
    Eigen::VectorXcd row_right_second;
};

struct NodePairs
{
    GlobalId g_el1, g_el2;
    LocalId l_el1, l_el2;
};

struct SBCNodeInfo
{
    NodePairs load;
    NodePairs unload;
};

using FieldComponentToFluxRows = std::map<std::pair<FieldType,Direction>,FluxRows>;
using ModalValues = Eigen::VectorXcd;
using NodalValues = Eigen::VectorXd;

class SBCSolver{
public:

    SBCSolver(Model&, ParFiniteElementSpace&, const SBCProperties*);

    void update(const Time& dt);
    void setTargetTime(const Time& t) { target_time = t; }
    void setPreTime(const Time& t) { pre_time = t; }
    void assignGlobalFields(Fields<ParFiniteElementSpace,ParGridFunction>* g_fields);

    
private:
    
    const SBCProperties* sbcp_;
    
    SBCNodeInfo dof_pair_;
    
    Model model_;
    
    Time target_time;
    Time pre_time = 0.0;
    Time dt_;
    
    SolverOptions opts_;

    FieldComponentToFluxRows nodal_to_modal_rows_;
    Eigen::MatrixXcd modal_to_nodal_matrix_;

    ModalValues modal_values_;
    NodalValues nodal_values_;

    std::vector<NodeId> target_ids_;
    
    Fields<ParFiniteElementSpace, ParGridFunction>* global_nodal_fields_;
    
    void findDoFPairs(Model&, ParFiniteElementSpace&);
    void initNodeIds(const std::pair<NodeId,NodeId>& el1_el2_ids);
    void loadNodalValuesAtFaces();
    void unloadNodalValuesAtFaces();
    void applyNodalToModalTransformation();
    void applyModalToNodalTransformation();

};

std::vector<NodeId> buildTargetNodeIds(size_t order, size_t num_of_segments);
Eigen::EigenSolver<Eigen::MatrixXd> applyEigenSolverOnGlobalOperator(const SparseMatrix& mat);
void updateModalValues(const FieldComponentToFluxRows& eigvecs, const Nodes& target_ids, const Fields<ParFiniteElementSpace, ParGridFunction>&, ModalValues& out);
void updateNodalValues(const Eigen::MatrixXcd& S, const ModalValues& q, NodalValues& x);
void loadEigenVectorFromOperator(const Eigen::MatrixXcd& op, const Nodes& target_ids, const size_t ndofs, FieldComponentToFluxRows& out);
ModalValues evolEigenvalueSystem(const ModalValues& q, const Eigen::VectorXcd& eigvals, const Nodes& target_ids, const Time dt);

}