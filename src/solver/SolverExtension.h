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

class SBCSolver{
public:

    SBCSolver(Model&, ParFiniteElementSpace&, const SBCProperties&);

    void setTargetTime(Time& t) { target_time = t; }
    void setPreTime(Time& t) { pre_time = t; }
    void assignGlobalFields(const Fields<ParFiniteElementSpace,ParGridFunction>* g_fields);
    NbrPairs getFieldValues(const FieldType f, const Direction d);
    
    private:
    
    SBCProperties sbcp_;
    
    std::unique_ptr<Mesh> mesh_;
    std::unique_ptr<ParMesh> pmesh_;
    std::unique_ptr<DG_FECollection> fec_;
    std::unique_ptr<ParFiniteElementSpace> fes_;
    std::vector<std::pair<NodeId, NodeId>> dof_pairs_;
    
    Model model_;
    
    Time target_time;
    Time pre_time = 0.0;
    Time dt_;
    
    SolverOptions opts_;
    
    const Fields<ParFiniteElementSpace, ParGridFunction>* global_fields_;
    Fields<ParFiniteElementSpace, ParGridFunction> sbc_fields_;
    
    std::unique_ptr<ODESolver> odeSolver_;
    std::unique_ptr<TimeDependentOperator> evolTDO_;
    
    void findDoFPairs(Model&, ParFiniteElementSpace&);
    
    void assignEvolutionOperator();

    void loadFieldValues(const FieldType, const Direction, const NbrPairs&);

};

class SBCTimeDependentOperator : public mfem::TimeDependentOperator
{
	public:
		static const int numberOfFieldComponents = 2;
		static const int numberOfMaxDimensions = 3;

		SBCTimeDependentOperator(Model&, ParFiniteElementSpace&);
		// virtual void Mult(const Vector& x, Vector& y) const;
		// void ImplicitSolve(const double dt, const Vector& x, Vector& k) override;

		const SparseMatrix& getConstGlobalOperator() { return *sbc_operator_.get(); }

	private:

		std::unique_ptr<mfem::SparseMatrix> sbc_operator_;

		ParFiniteElementSpace& fes_;
		Model& model_;

		mutable std::array<ParGridFunction, 3> eOld_, hOld_;
};

}