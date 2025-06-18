#pragma once

#include "ProblemDescription.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolutionMethods.h"
#include "evolution/MaxwellEvolutionMethods.h"

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"

#include "Types.h"

#include <chrono>
#include <iostream>
#include <string>

namespace maxwell
{

	using namespace mfem;
	using namespace mfemExtension;

	static FluxBdrCoefficientsCentered bdrCentCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
	};

	static FluxBdrCoefficientsUpwind bdrUpwindCoeff{
		{BdrCond::PEC, {2.0, 0.0}},
		{BdrCond::PMC, {0.0, 2.0}},
		{BdrCond::SMA, {1.0, 1.0}},
		{BdrCond::SurfaceCond, {1.0, 1.0}},
	};

	static FluxSrcCoefficientsCentered srcCentCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
	};

	static FluxSrcCoefficientsUpwind srcUpwindCoeff{
		{BdrCond::TotalFieldIn, {1.0, 1.0}},
	};

	static FieldType altField(const FieldType &f)
	{
		switch (f)
		{
		case FieldType::E:
			return FieldType::H;
		case FieldType::H:
			return FieldType::E;
		default:
			throw std::runtime_error("Incorrect FieldType in input.");
		}
	}

	static void loadBlockInGlobalAtIndices(const SparseMatrix &blk, SparseMatrix &dst, const std::pair<Array<int>, Array<int>> &ids, const double fieldSign)
	{
		MFEM_ASSERT(blk.NumRows() == ids.first.Size(), "Block Sparse NumRows does not match intended number of Rows.");
		MFEM_ASSERT(blk.NumCols() == ids.second.Size(), "Block Sparse NumCols does not match intended number of Cols.");
		Array<int> cols;
		Vector vals;
		for (auto r{0}; r < ids.first.Size(); r++)
		{
			blk.GetRow(r, cols, vals);
			for (auto c{0}; c < cols.Size(); c++)
			{
				dst.Add(ids.first[r], ids.second[cols[c]], vals[c] * fieldSign);
			}
		}
	}

	struct GlobalIndices
	{
		GlobalIndices(const int blockSize)
		{
			for (auto f : {E, H})
			{
				for (auto d : {X, Y, Z})
				{
					index[f][d] = Array<int>(blockSize);
					for (auto i{0}; i < blockSize; i++)
					{
						index[f][d][i] = i + (3 * f + d) * blockSize;
					}
				}
			}
		}

		std::array<std::array<mfem::Array<int>, 3>, 2> index;
	};

	static std::map<BdrCond, std::vector<double>> bdrCoeffCheck(double alpha)
	{
		std::map<BdrCond, std::vector<double>> res;
		if (alpha == 0.0)
		{
			res = bdrCentCoeff;
		}
		else
		{
			res = bdrUpwindCoeff;
		}
		return res;
	}

	template <typename FES, typename BF>
	std::unique_ptr<BF> buildByMult(
		const SparseMatrix &op1,
		const SparseMatrix &op2,
		FES &fes)
	{
		SparseMatrix *matrix = mfem::Mult(op1, op2);

		std::unique_ptr<BF> res = std::make_unique<BF>(&fes);

		res->Assemble();
		res->Finalize();
		res->SpMat().Swap(*matrix);
		delete matrix;

		return res;
	}

	void loadBlockInGlobalAtIndices(const SparseMatrix &blk, SparseMatrix &dst, const std::pair<Array<int>, Array<int>> &ids, const double fieldSign = 1.0);

	template <typename FES>
	class DGOperatorFactory
	{
	public:
		DGOperatorFactory(ProblemDescription &pd, FES &fes);

		// Methods for speficic FieldType or Direction Operators //
		template <typename BF>
		std::unique_ptr<BF> buildInverseMassMatrixSubOperator(const FieldType &f);

		template <typename BF>
		std::unique_ptr<BF> buildDerivativeSubOperator(const Direction &d);
		template <typename BF>
		std::unique_ptr<BF> buildZeroNormalSubOperator(const FieldType &f);
		template <typename BF>
		std::unique_ptr<BF> buildOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);
		template <typename BF>
		std::unique_ptr<BF> buildTwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);

		template <typename BF>
		std::unique_ptr<BF> buildZeroNormalIBFISubOperator(const FieldType &f);
		template <typename BF>
		std::unique_ptr<BF> buildOneNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);
		template <typename BF>
		std::unique_ptr<BF> buildTwoNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms);

		// Methods for complete Maxwell Operators //
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellInverseMassMatrixOperator();
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellTFSFInverseMassMatrixOperator();

		template <typename BF>
		std::array<std::array<std::unique_ptr<BF>, 3>, 2> buildMaxwellDirectionalOperator();
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellZeroNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> buildMaxwellOneNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> buildMaxwellTwoNormalOperator();

		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildMaxwellIntBdrZeroNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> buildMaxwellIntBdrOneNormalOperator();
		template <typename BF>
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> buildMaxwellIntBdrTwoNormalOperator();

		// Methors for complete Global Operators //
		template <typename BF>
		std::array<std::unique_ptr<BF>, 2> buildGlobalInverseMassMatrixOperator();

		// void addGlobalZeroNormalIBFIOperators(SparseMatrix* global);
		// void addGlobalOneNormalIBFIOperators(SparseMatrix* global);
		// void addGlobalTwoNormalIBFIOperators(SparseMatrix* global);
		// void addGlobalDirectionalOperators(SparseMatrix* global);
		// void addGlobalZeroNormalOperators(SparseMatrix* global);
		// void addGlobalOneNormalOperators(SparseMatrix* global);
		// void addGlobalTwoNormalOperators(SparseMatrix* global);

		std::unique_ptr<SparseMatrix> buildTFSFGlobalOperator();
		std::unique_ptr<SparseMatrix> buildGlobalOperator();

	private:
		ProblemDescription pd_;
		FES fes_;
	};

	template <typename FES>
	DGOperatorFactory<FES>::DGOperatorFactory(ProblemDescription &pd, FES &fes) : pd_(pd),
																				  fes_(fes)
	{
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildInverseMassMatrixSubOperator(const FieldType &f)
	{
		Vector aux{pd_.model.buildEpsMuPiecewiseVector(f)};
		PWConstCoefficient PWCoeff(aux);

		auto res = std::make_unique<BF>(&fes_);
		res->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildDerivativeSubOperator(const Direction &d)
	{
		auto res = std::make_unique<BF>(&fes_);

		if (d >= fes_.GetMesh()->Dimension())
		{
			res->Assemble();
			res->Finalize();
			return res;
		}

		ConstantCoefficient coeff = (d <= fes_.GetMesh()->Dimension()) ? ConstantCoefficient(1.0) : ConstantCoefficient(0.0);
		res->AddDomainIntegrator(
			new DerivativeIntegrator(coeff, d));

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildZeroNormalSubOperator(const FieldType &f)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildOneNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildTwoNormalSubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size())
		{
			for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
			{
				res->AddInteriorFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), kv.second);
			}
		}
		else
		{
			res->AddInteriorFaceIntegrator(
				new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha));
		}

		for (auto &kv : pd_.model.getBoundaryToMarker())
		{

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA)
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else
			{
				res->AddBdrFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildZeroNormalIBFISubOperator(const FieldType &f)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator({}, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildOneNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::unique_ptr<BF> DGOperatorFactory<FES>::buildTwoNormalIBFISubOperator(const FieldType &f, const std::vector<Direction> &dirTerms)
	{
		auto res = std::make_unique<BF>(&fes_);

		for (auto &kv : pd_.model.getInteriorBoundaryToMarker())
		{
			if (kv.first != BdrCond::TotalFieldIn)
			{
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first)
				{
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, 1.0), kv.second);
					break;
				default:
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
					break;
				}
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellInverseMassMatrixOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			res[f] = buildInverseMassMatrixSubOperator<BF>(f);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::unique_ptr<BF>, 3>, 2> DGOperatorFactory<FES>::buildMaxwellDirectionalOperator()
	{
		std::array<std::array<std::unique_ptr<BF>, 3>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				res[f][d] = buildByMult<FES, BF>(MInv->SpMat(), buildDerivativeSubOperator<BF>(d)->SpMat(), fes_);
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellZeroNormalOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			res[f] = buildByMult<FES, BF>(MInv->SpMat(), buildZeroNormalSubOperator<BF>(f)->SpMat(), fes_);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellOneNormalOperator()
	{
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					res[f][f2][d] = buildByMult<FES, BF>(MInv->SpMat(), buildOneNormalSubOperator<BF>(f2, {d})->SpMat(), fes_);
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					for (auto d2{X}; d2 <= Z; d2++)
					{
						res[f][f2][d][d2] = buildByMult<FES, BF>(MInv->SpMat(), buildTwoNormalSubOperator<BF>(f2, {d, d2})->SpMat(), fes_);
					}
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrZeroNormalOperator()
	{
		std::array<std::unique_ptr<BF>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			res[f] = buildByMult<FES, BF>(MInv->SpMat(), buildZeroNormalIBFISubOperator<BF>(f)->SpMat(), fes_);
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrOneNormalOperator()
	{
		std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					res[f][f2][d] = buildByMult<FES, BF>(MInv->SpMat(), buildOneNormalIBFISubOperator<BF>(f2, {d})->SpMat(), fes_);
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> DGOperatorFactory<FES>::buildMaxwellIntBdrTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<std::unique_ptr<BF>, 3>, 3>, 2>, 2> res;
		for (auto f : {E, H})
		{
			auto MInv = buildInverseMassMatrixSubOperator<BF>(f);
			for (auto d{X}; d <= Z; d++)
			{
				for (auto f2 : {E, H})
				{
					for (auto d2{X}; d2 <= Z; d2++)
					{
						res[f][f2][d][d2] = buildByMult<FES, BF>(MInv->SpMat(), buildTwoNormalIBFISubOperator<BF>(f2, {d, d2})->SpMat(), fes_);
					}
				}
			}
		}
		return res;
	}

	template <typename FES>
	template <typename BF>
	std::array<std::unique_ptr<BF>, 2> DGOperatorFactory<FES>::buildGlobalInverseMassMatrixOperator()
	{
		return buildMaxwellInverseMassMatrixOperator();
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildTFSFGlobalOperator()
	{

		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * fes_.GetNDofs());

		res->Finalize();
		return res;
	}

	template <typename FES>
	std::unique_ptr<SparseMatrix> DGOperatorFactory<FES>::buildGlobalOperator()
	{

		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * fes_.GetNDofs());

		res->Finalize();
		return res;
	}

}