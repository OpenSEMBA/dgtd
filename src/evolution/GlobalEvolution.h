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

		const SparseMatrix& getConstGlobalOperator() { return *globalOperator_.get(); }

	private:

		std::unique_ptr<mfem::SparseMatrix> globalOperator_;
		std::unique_ptr<mfem::SparseMatrix> TFSFOperator_;
		Array<int> sub_to_parent_ids_;

		mfem::ParFiniteElementSpace& fes_;
		Model& model_;
		SourcesManager& srcmngr_;
		EvolutionOptions& opts_;

		mutable std::array<ParGridFunction, 3> eOld_, hOld_;

	};

mfem::Vector load_tfsf_into_single_vector_gpu(const FieldGridFuncs& func);

void load_tfsf_into_out_vector_gpu(const mfem::Array<int>& sub_to_parent_ids_, 
                                   const mfem::Vector& tempTFSF,                               
                                         mfem::Vector& out,
                                   const int global_ndofs,
                                   const int tfsf_ndofs);

FieldGridFuncs eval_time_var_field_gpu(const Time time, SourcesManager& sm);

}