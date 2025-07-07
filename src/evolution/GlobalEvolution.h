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

		mutable std::array<mfem::ParGridFunction, 3> eOld_, hOld_, eFull_, hFull_;
		mutable mfem::Vector inNew_;

	};

void load_in_to_eh_gpu(mfem::Vector& in, 
                        std::array<mfem::Vector, 3>& eOld, 
                        std::array<mfem::Vector, 3>& hOld, 
                        const int ndofs);

void load_eh_to_innew_gpu(const std::array<mfem::Vector, 3>& eOld,
                          const std::array<mfem::Vector, 3>& hOld,
                          mfem::Vector& inNew,
                          int ndofs,
                          int nbrSize);

void load_nbr_to_innew_gpu(const std::array<mfem::ParGridFunction, 3>& eOldNbr,
                   const std::array<mfem::ParGridFunction, 3>& hOldNbr,
                   mfem::Vector& inNew,
                   const int ndofs,
                   const int nbrSize);

}