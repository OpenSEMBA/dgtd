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

    bool implicit_ode = false;

    SBCProperties(size_t segnum, size_t o, size_t mat_w, bool is_implicit = false) : 
    num_of_segments(segnum), order(o), material_width(mat_w), implicit_ode(is_implicit){}

};

class SBCSolver{
public:

    SBCSolver(Model&, FiniteElementSpace&, const SBCProperties&);

    void setGlobalTime(Time& t) { global_time_ = t; }
    void setSBCTime(Time& t) { sbc_time_ = t; }

private:

    SBCProperties sbcp_;

    std::unique_ptr<Mesh> mesh_;
    std::unique_ptr<DG_FECollection> fec_;
    std::unique_ptr<FiniteElementSpace> fes_;
    std::vector<std::pair<NodeId, NodeId>> dof_pairs_;

    Time global_time_;
    Time sbc_time_;
    Time dt_;

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