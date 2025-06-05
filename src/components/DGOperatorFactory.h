#pragma once

#include "ProblemDefinition.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "evolution/MaxwellEvolutionMethods.h"

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

#include "Types.h"

namespace maxwell {

	using FiniteElementOperator = std::unique_ptr<mfemExtension::BilinearForm>;

	FieldType altField(const FieldType& f);

	FiniteElementOperator buildByMult(const SparseMatrix& op1, const SparseMatrix& op2, FiniteElementSpace& fes);

	struct GlobalIndices {
		GlobalIndices(const int blockSize);
		std::array<std::array<mfem::Array<int>, 3>, 2> index;
	};

	void loadBlockInGlobalAtIndices(const SparseMatrix& blk, SparseMatrix& dst, const std::pair<Array<int>, Array<int>>& ids, const double fieldSign = 1.0, bool temp_dbg = false);

	class DGOperatorFactory {
	public:

		DGOperatorFactory(ProblemDescription& pd, mfem::FiniteElementSpace& fes);

		// Methods for speficic FieldType or Direction Operators //

		FiniteElementOperator buildInverseMassMatrixSubOperator(const FieldType& f);
		FiniteElementOperator buildDerivativeSubOperator(const Direction& d);
		FiniteElementOperator buildZeroNormalSubOperator(const FieldType& f);
		FiniteElementOperator buildOneNormalSubOperator(const FieldType& f, const std::vector<Direction>& dirTerms);
		FiniteElementOperator buildTwoNormalSubOperator(const FieldType& f, const std::vector<Direction>& dirTerms);

		FiniteElementOperator buildZeroNormalIBFISubOperator(const FieldType& f);
		FiniteElementOperator buildOneNormalIBFISubOperator(const FieldType& f, const std::vector<Direction>& dirTerms);
		FiniteElementOperator buildTwoNormalIBFISubOperator(const FieldType& f, const std::vector<Direction>& dirTerms);

		// Methods for complete Maxwell Operators //

		std::array<FiniteElementOperator, 2> buildMaxwellInverseMassMatrixOperator();
		std::array<FiniteElementOperator, 2> buildMaxwellTFSFInverseMassMatrixOperator();

		std::array<std::array<FiniteElementOperator, 3>, 2> buildMaxwellDirectionalOperator();
		std::array<FiniteElementOperator, 2> buildMaxwellZeroNormalOperator();
		std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> buildMaxwellOneNormalOperator();
		std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> buildMaxwellTwoNormalOperator();

		std::array<FiniteElementOperator, 2> buildMaxwellIntBdrZeroNormalOperator();
		std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> buildMaxwellIntBdrOneNormalOperator();
		std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> buildMaxwellIntBdrTwoNormalOperator();

		// Methors for complete Global Operators //

		std::array<FiniteElementOperator, 2> buildGlobalInverseMassMatrixOperator();
		void addGlobalZeroNormalIBFIOperators(SparseMatrix* global);
		void addGlobalOneNormalIBFIOperators(SparseMatrix* global);
		void addGlobalTwoNormalIBFIOperators(SparseMatrix* global);
		void addGlobalDirectionalOperators(SparseMatrix* global);
		void addGlobalZeroNormalOperators(SparseMatrix* global);
		void addGlobalOneNormalOperators(SparseMatrix* global);
		void addGlobalTwoNormalOperators(SparseMatrix* global);

		std::unique_ptr<SparseMatrix> buildTFSFGlobalOperator();
		std::unique_ptr<SparseMatrix> buildGlobalOperator();

	private:

		ProblemDescription pd_;
		mfem::FiniteElementSpace fes_;
		bool temp_dbg = false;
};

}