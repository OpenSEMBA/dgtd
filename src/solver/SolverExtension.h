#include "components/Model.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "SolverOptions.h"

#include <mfem.hpp>

namespace maxwell 
{

using namespace mfem;

class SBCSolver{
public:

    SBCSolver(Model&, FiniteElementSpace&, const SBCProperties&);

    void setGlobalTime(Time& t) { global_time_ = t; }
    void setSBCTime(Time& t) { sbc_time_ = t; }
    void resetFields();
    void assignGlobalFields(const Fields<ParFiniteElementSpace,ParGridFunction>* g_fields);
    std::pair<double, double> getFieldPairAfterCalculation(const FieldType f, const Direction d);
    
    private:
    
    SBCProperties sbcp_;
    
    std::unique_ptr<Mesh> mesh_;
    std::unique_ptr<DG_FECollection> fec_;
    std::unique_ptr<FiniteElementSpace> fes_;
    std::vector<std::pair<NodeId, NodeId>> dof_pairs_;

    Time global_time_;
    Time sbc_time_;
    Time dt_;

    const Fields<ParFiniteElementSpace,ParGridFunction>* global_fields_;
    Fields<FiniteElementSpace,GridFunction> sbc_fields_;
    
    std::unique_ptr<ODESolver> odeSolver_;
    std::unique_ptr<TimeDependentOperator> evolTDO_;
    
    void findDoFPairs(Model&, FiniteElementSpace&);
    
    void assignODESolver();
    void assignEvolutionOperator();

    void estimateTimeStep();

};

class SBCTimeDependentOperator : public mfem::TimeDependentOperator
{
	public:
		static const int numberOfFieldComponents = 2;
		static const int numberOfMaxDimensions = 3;

		SBCTimeDependentOperator(Model&, FiniteElementSpace&);
		virtual void Mult(const Vector& x, Vector& y) const;
		void ImplicitSolve(const double dt, const Vector& x, Vector& k) override;

		const SparseMatrix& getConstGlobalOperator() { return *sbc_operator_.get(); }

	private:

		std::unique_ptr<mfem::SparseMatrix> sbc_operator_;

		FiniteElementSpace& fes_;
		Model& model_;

		mutable std::array<ParGridFunction, 3> eOld_, hOld_;
};

}