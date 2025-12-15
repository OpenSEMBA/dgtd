#pragma once

#include "solver/SourcesManager.h"
#include "components/DGOperatorFactory.h"
#include "components/Types.h"

namespace maxwell {

	class GlobalEvolution : public mfem::TimeDependentOperator {
	public:
		static const int numberOfFieldComponents = 2;
		static const int numberOfMaxDimensions = 3;

		GlobalEvolution(mfem::ParFiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
		virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;
		void ImplicitSolve(const double dt, const mfem::Vector& x, mfem::Vector& k) override;

		const SparseMatrix& getConstGlobalOperator() { return *globalOperator_.get(); }

	private:

		std::unique_ptr<mfem::SparseMatrix> globalOperator_;
		std::unique_ptr<mfem::SparseMatrix> TFSFOperator_;
		std::unique_ptr<mfem::SparseMatrix> SGBCOperator_;
		Array<int> tfsf_sub_to_parent_ids_;
		Array<int> sgbc_sub_to_parent_ids_;

		mfem::ParFiniteElementSpace& fes_;
		Model& model_;
		SourcesManager& srcmngr_;
		EvolutionOptions& opts_;

		mutable std::array<ParGridFunction, 3> eOld_, hOld_;

	};


void load_in_to_eh_gpu(const mfem::Vector& in, 
                        std::array<ParGridFunction, 3>& e,
                        std::array<ParGridFunction, 3>& h,
                        int ndofs);
	
void load_eh_to_innew_gpu(const mfem::Vector& inOld,
                          mfem::Vector& inNew,
                          int ndofs,
                          int nbrSize);

void load_nbr_to_innew_gpu(const std::array<mfem::ParGridFunction, 3>& eOldNbr,
                   const std::array<mfem::ParGridFunction, 3>& hOldNbr,
                   mfem::Vector& inNew,
                   const int ndofs,
                   const int nbrSize);

mfem::Vector load_tfsf_into_single_vector_gpu(const FieldGridFuncs& func);

void load_tfsf_into_out_vector_gpu(const mfem::Array<int>& tfsf_sub_to_parent_ids_, 
                                   const mfem::Vector& tempTFSF,                               
                                         mfem::Vector& out,
                                   const int global_ndofs,
                                   const int tfsf_ndofs);

FieldGridFuncs eval_time_var_field_gpu(const Time time, SourcesManager& sm);

}