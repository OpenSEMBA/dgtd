#pragma once

#include "solver/SourcesManager.h"
#include "components/DGOperatorFactory.h"
#include "components/Types.h"

namespace maxwell {

	class GlobalEvolution : public mfem::TimeDependentOperator {
	public:
		static const int numberOfFieldComponents = 2;
		static const int numberOfMaxDimensions = 3;

		GlobalEvolution(mfem::FiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
		virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

		const SparseMatrix& getConstGlobalOperator() { return *globalOperator_.get(); }

	private:

		std::unique_ptr<mfem::SparseMatrix> globalOperator_;
		std::unique_ptr<mfem::SparseMatrix> TFSFOperator_;
		Array<int> sub_to_parent_ids_;

		mfem::FiniteElementSpace& fes_;
		Model& model_;
		SourcesManager& srcmngr_;
		EvolutionOptions& opts_;

	};

}