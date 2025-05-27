#include "DGOperatorFactory.h"

#include <chrono>
#include <iostream>
#include <string>

namespace maxwell {
	
	using namespace mfem;
	using namespace mfemExtension;

	FluxBdrCoefficientsCentered bdrCentCoeff{
		{BdrCond::PEC               ,	{2.0, 0.0}},
		{BdrCond::PMC               ,	{0.0, 2.0}},
		{BdrCond::SMA               ,	{1.0, 1.0}},
		{BdrCond::SurfaceCond       ,   {1.0, 1.0}},
	};

	FluxBdrCoefficientsUpwind bdrUpwindCoeff{
		{BdrCond::PEC               ,	{2.0, 0.0}},
		{BdrCond::PMC               ,	{0.0, 2.0}},
		{BdrCond::SMA               ,	{1.0, 1.0}},
		{BdrCond::SurfaceCond       ,   {1.0, 1.0}},
	};

	FluxSrcCoefficientsCentered srcCentCoeff{
		{BdrCond::TotalFieldIn      ,	{1.0, 1.0}},
	};

	FluxSrcCoefficientsUpwind srcUpwindCoeff{
		{BdrCond::TotalFieldIn      ,	{1.0, 1.0}},
	};

	FieldType altField(const FieldType& f)
	{
		switch (f) {
		case FieldType::E:
			return FieldType::H;
		case FieldType::H:
			return FieldType::E;
		default:
			throw std::runtime_error("Incorrect FieldType in input.");
		}
	}

	std::map<BdrCond, std::vector<double>> bdrCoeffCheck(double alpha)
	{
		std::map<BdrCond, std::vector<double>> res;
		if (alpha == 0.0) {
			res = bdrCentCoeff;
		}
		else{
			res = bdrUpwindCoeff;
		}
		return res;
	}

	FiniteElementOperator buildByMult(
		const BilinearForm& op1,
		const BilinearForm& op2,
		FiniteElementSpace& fes)
	{
		auto aux = std::unique_ptr<SparseMatrix>(mfem::Mult(op1.SpMat(), op2.SpMat()));
		auto res = std::make_unique<BilinearForm>(&fes);
		res->Assemble();
		res->Finalize();
		res->SpMat().Swap(*aux);

		return res;
	}

	void loadBlockInGlobalAtIndices(const SparseMatrix& blk, SparseMatrix& dst, const std::pair<Array<int>, Array<int>>& ids, const double fieldSign)
	{
		MFEM_ASSERT(blk.NumRows() == ids.first.Size(),  "Block Sparse NumRows does not match intended number of Rows.");
		MFEM_ASSERT(blk.NumCols() == ids.second.Size(), "Block Sparse NumCols does not match intended number of Cols.");
		Array<int> cols;
		Vector vals;
		for (auto r{ 0 }; r < ids.first.Size(); r++) {
			blk.GetRow(r, cols, vals);
			for (auto c{ 0 }; c < cols.Size(); c++) {
				dst.Add(ids.first[r], ids.second[cols[c]], vals[c] * fieldSign);
			}
		}
	}

	GlobalIndices::GlobalIndices(const int blockSize)
	{
		for (auto f : { E, H }) {
			for (auto d : { X, Y, Z }) {
				index[f][d] = Array<int>(blockSize);
				for (auto i{ 0 }; i < blockSize; i++) {
					index[f][d][i] = i + (3 * f + d) * blockSize;
				}
			}
		}
	}

	DGOperatorFactory::DGOperatorFactory(ProblemDescription& pd, FiniteElementSpace& fes) :
		pd_(pd),
		fes_(fes)
	{}

	FiniteElementOperator DGOperatorFactory::buildInverseMassMatrixSubOperator(const FieldType& f)
	{
		Vector aux{ pd_.model.buildEpsMuPiecewiseVector(f) };
		PWConstCoefficient PWCoeff(aux);

		auto res = std::make_unique<BilinearForm>(&fes_);
		res->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));

		res->Assemble();
		res->Finalize();
		return res;
	}

	FiniteElementOperator DGOperatorFactory::buildDerivativeSubOperator(const Direction& d)
	{
		auto res = std::make_unique<BilinearForm>(&fes_);

		if (d >= fes_.GetMesh()->Dimension()) {
			res->Assemble();
			res->Finalize();
			return res;
		}

		ConstantCoefficient coeff = (d <= fes_.GetMesh()->Dimension()) ? ConstantCoefficient(1.0) : ConstantCoefficient(0.0);
		res->AddDomainIntegrator(
			new DerivativeIntegrator(coeff, d)
		);

		res->Assemble();
		res->Finalize();
		return res;
	}

	FiniteElementOperator DGOperatorFactory::buildZeroNormalSubOperator(const FieldType& f)
	{
		auto res = std::make_unique<BilinearForm>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size()) {
			for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
				res->AddInteriorFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha), kv.second);
			}
		}
		else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGZeroNormalJumpIntegrator(pd_.opts.alpha));
		}

		for (auto& kv : pd_.model.getBoundaryToMarker()) {

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA) {
				res->AddBdrFaceIntegrator(
					new MaxwellDGZeroNormalJumpIntegrator(c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else {
				res->AddBdrFaceIntegrator(
					new mfemExtension::MaxwellSMAJumpIntegrator({}), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}

	FiniteElementOperator DGOperatorFactory::buildOneNormalSubOperator(const FieldType& f, const std::vector<Direction>& dirTerms)
	{
		auto res = std::make_unique<BilinearForm>(&fes_);
		ConstantCoefficient one(1.0);
		if (pd_.model.getInteriorBoundaryToMarker().size()) {
			for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
				res->AddInteriorFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0), kv.second);
			}
		}
		else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGOneNormalJumpIntegrator(dirTerms, 1.0));
		}

		for (auto& kv : pd_.model.getBoundaryToMarker()) {

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA) {
				res->AddBdrFaceIntegrator(
					new MaxwellDGOneNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			}
			else {
				res->AddBdrFaceIntegrator(
					new mfemExtension::MaxwellSMAJumpIntegrator(dirTerms), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}
	FiniteElementOperator DGOperatorFactory::buildTwoNormalSubOperator(const FieldType& f, const std::vector<Direction>& dirTerms)
	{
		auto res = std::make_unique<BilinearForm>(&fes_);
		if (pd_.model.getInteriorBoundaryToMarker().size()) {
			for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
				res->AddInteriorFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha), kv.second);
			}
		}
		else {
			res->AddInteriorFaceIntegrator(
				new MaxwellDGTwoNormalJumpIntegrator(dirTerms, pd_.opts.alpha));
		}

		for (auto& kv : pd_.model.getBoundaryToMarker()) {

			auto c = bdrCoeffCheck(pd_.opts.alpha);
			if (kv.first != BdrCond::SMA) {
				res->AddBdrFaceIntegrator(
					new MaxwellDGTwoNormalJumpIntegrator(dirTerms, c[kv.first].at(f) * pd_.opts.alpha), kv.second);
			}
			else {
				res->AddBdrFaceIntegrator(
					new mfemExtension::MaxwellSMAJumpIntegrator(dirTerms), kv.second);
			}
		}

		res->Assemble();
		res->Finalize();
		return res;
	}
	FiniteElementOperator DGOperatorFactory::buildZeroNormalIBFISubOperator(const FieldType& f)
	{
		auto res = std::make_unique<mfemExtension::BilinearForm>(&fes_);

		for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
			if (kv.first != BdrCond::TotalFieldIn) {
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first) {
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellSMAJumpIntegrator({}), kv.second);
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
	FiniteElementOperator DGOperatorFactory::buildOneNormalIBFISubOperator(const FieldType& f, const std::vector<Direction>& dirTerms)
	{
		auto res = std::make_unique<mfemExtension::BilinearForm>(&fes_);

		for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
			if (kv.first != BdrCond::TotalFieldIn) {
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first) {
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellSMAJumpIntegrator(dirTerms), kv.second);
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
	FiniteElementOperator DGOperatorFactory::buildTwoNormalIBFISubOperator(const FieldType& f, const std::vector<Direction>& dirTerms)
	{
		auto res = std::make_unique<mfemExtension::BilinearForm>(&fes_);

		for (auto& kv : pd_.model.getInteriorBoundaryToMarker()) {
			if (kv.first != BdrCond::TotalFieldIn) {
				auto c = bdrCoeffCheck(pd_.opts.alpha);
				switch (kv.first) {
				case (BdrCond::SMA):
					res->AddInternalBoundaryFaceIntegrator(
						new mfemExtension::MaxwellSMAJumpIntegrator(dirTerms), kv.second);
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

	std::array<FiniteElementOperator, 2> DGOperatorFactory::buildMaxwellInverseMassMatrixOperator()
	{
		std::array<FiniteElementOperator, 2> res;
		for (auto f : { E, H }) {
			res[f] = buildInverseMassMatrixSubOperator(f);
		}
		return res;
	}
	std::array<std::array<FiniteElementOperator, 3>, 2> DGOperatorFactory::buildMaxwellDirectionalOperator()
	{
		std::array<std::array<FiniteElementOperator, 3>, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			for (auto d{ X }; d <= Z; d++) {
				res[f][d] = buildByMult(*MInv, *buildDerivativeSubOperator(d), fes_);
			}
		}
		return res;
	}
	std::array<FiniteElementOperator, 2> DGOperatorFactory::buildMaxwellZeroNormalOperator()
	{
		std::array<FiniteElementOperator, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			res[f] = buildByMult(*MInv, *buildZeroNormalSubOperator(f), fes_);
		}
		return res;
	}
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> DGOperatorFactory::buildMaxwellOneNormalOperator()
	{
		std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					res[f][f2][d] = buildByMult(*MInv, *buildOneNormalSubOperator(f2, { d }), fes_);
				}
			}
		}
		return res;
	}
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> DGOperatorFactory::buildMaxwellTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						res[f][f2][d][d2] = buildByMult(*MInv, *buildTwoNormalSubOperator(f2, { d, d2 }), fes_);
					}
				}
			}
		}
		return res;
	}
	std::array<FiniteElementOperator, 2> DGOperatorFactory::buildMaxwellIntBdrZeroNormalOperator()
	{
		std::array<FiniteElementOperator, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			res[f] = buildByMult(*MInv, *buildZeroNormalIBFISubOperator(f), fes_);
		}
		return res;
	}
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> DGOperatorFactory::buildMaxwellIntBdrOneNormalOperator()
	{
		std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					res[f][f2][d] = buildByMult(*MInv, *buildOneNormalIBFISubOperator(f2, { d }), fes_);
				}
			}
		}
		return res;
	}
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> DGOperatorFactory::buildMaxwellIntBdrTwoNormalOperator()
	{
		std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> res;
		for (auto f : { E, H }) {
			auto MInv = buildInverseMassMatrixSubOperator(f);
			for (auto d{ X }; d <= Z; d++) {
				for (auto f2 : { E, H }) {
					for (auto d2{ X }; d2 <= Z; d2++) {
						res[f][f2][d][d2] = buildByMult(*MInv, *buildTwoNormalIBFISubOperator(f2, { d, d2 }), fes_);
					}
				}
			}
		}
		return res;
	}


	std::array<FiniteElementOperator, 2> DGOperatorFactory::buildGlobalInverseMassMatrixOperator()
	{
		return buildMaxwellInverseMassMatrixOperator();
	}
	void DGOperatorFactory::addGlobalZeroNormalIBFIOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			auto op = buildByMult(*MInv[f], *buildZeroNormalIBFISubOperator(f), fes_);
			for (auto d : { X, Y, Z }) {
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][d], globalId.index[f][d]),
					-1.0
				);
			}
		}
	}
	void DGOperatorFactory::addGlobalOneNormalIBFIOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult(*MInv[f], *buildOneNormalIBFISubOperator(altField(f), { x }), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][y], globalId.index[altField(f)][z]),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][z], globalId.index[altField(f)][y]),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}
	void DGOperatorFactory::addGlobalTwoNormalIBFIOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					auto op = buildByMult(*MInv[f], *buildTwoNormalIBFISubOperator(f, { d, d2 }), fes_);
					loadBlockInGlobalAtIndices(
						op->SpMat(),
						*global,
						std::make_pair(globalId.index[f][d], globalId.index[f][d2])
					);
				}
			}
		}
	}
	void DGOperatorFactory::addGlobalDirectionalOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult(*MInv[f], *buildDerivativeSubOperator(x), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][z], globalId.index[altField(f)][y]),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][y], globalId.index[altField(f)][z]),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}
	void DGOperatorFactory::addGlobalZeroNormalOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			auto op = buildByMult(*MInv[f], *buildZeroNormalSubOperator(f), fes_);
			for (auto d : { X, Y, Z }) {
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][d], globalId.index[f][d]),
					-1.0
				);
			}
		}
	}
	void DGOperatorFactory::addGlobalOneNormalOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			for (auto x{ X }; x <= Z; x++) {
				auto y = (x + 1) % 3;
				auto z = (x + 2) % 3;
				auto op = buildByMult(*MInv[f], *buildOneNormalSubOperator(altField(f), { x }), fes_);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][y], globalId.index[altField(f)][z]),
					1.0 - double(f) * 2.0
				);
				loadBlockInGlobalAtIndices(
					op->SpMat(),
					*global,
					std::make_pair(globalId.index[f][z], globalId.index[altField(f)][y]),
					-1.0 + double(f) * 2.0
				);
			}
		}
	}
	void DGOperatorFactory::addGlobalTwoNormalOperators(SparseMatrix* global)
	{
		GlobalIndices globalId(fes_.GetNDofs());
		for (auto f : { E, H }) {
			auto MInv = buildGlobalInverseMassMatrixOperator();
			for (auto d{ X }; d <= Z; d++) {
				for (auto d2{ X }; d2 <= Z; d2++) {
					auto op = buildByMult(*MInv[f], *buildTwoNormalSubOperator(f, { d, d2 }), fes_);
					loadBlockInGlobalAtIndices(
						op->SpMat(),
						*global,
						std::make_pair(globalId.index[f][d], globalId.index[f][d2])
					);
				}
			}
		}
	}

	std::unique_ptr<SparseMatrix> DGOperatorFactory::buildTFSFGlobalOperator()
	{

#ifdef SHOW_TIMER_INFORMATION
		auto startTime{ std::chrono::high_resolution_clock::now() };
#endif

		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * fes_.GetNDofs());


#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling TFSF Inverse Mass One-Normal Operators" << std::endl;
#endif

		addGlobalOneNormalOperators(res.get());

		if (pd_.opts.alpha > 0.0) {

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Zero-Normal Operators" << std::endl;
#endif	

			addGlobalZeroNormalOperators(res.get());


#ifdef SHOW_TIMER_INFORMATION	
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling TFSF Inverse Mass Two-Normal Operators" << std::endl;
#endif	

			addGlobalTwoNormalOperators(res.get());

		}
		
		res->Finalize();
		return res;

	}

	std::unique_ptr<SparseMatrix> DGOperatorFactory::buildGlobalOperator()
	{

#ifdef SHOW_TIMER_INFORMATION
		auto startTime{ std::chrono::high_resolution_clock::now() };
#endif
		std::unique_ptr<SparseMatrix> res = std::make_unique<SparseMatrix>(6 * fes_.GetNDofs(), 6 * fes_.GetNDofs());

		if (pd_.model.getInteriorBoundaryToMarker().size() != 0) { //IntBdrConds

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling IBFI Inverse Mass One-Normal Operators" << std::endl;
#endif

			addGlobalOneNormalIBFIOperators(res.get());

			if (pd_.opts.alpha != 0.0) {

#ifdef SHOW_TIMER_INFORMATION
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				std::cout << "Assembling IBFI Inverse Mass Zero-Normal Operators" << std::endl;
#endif

				addGlobalZeroNormalIBFIOperators(res.get());

#ifdef SHOW_TIMER_INFORMATION
				std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
				std::cout << "Assembling IBFI Inverse Mass Two-Normal Operators" << std::endl;
#endif

				addGlobalTwoNormalIBFIOperators(res.get());

			}
		}

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass Stiffness Operators" << std::endl;
#endif

		addGlobalDirectionalOperators(res.get());

#ifdef SHOW_TIMER_INFORMATION
		std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
		std::cout << "Assembling Standard Inverse Mass One-Normal Operators" << std::endl;
#endif

		addGlobalOneNormalOperators(res.get());

		if (pd_.opts.alpha != 0.0) {

#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Zero-Normal Operators" << std::endl;
#endif
			addGlobalZeroNormalOperators(res.get());


#ifdef SHOW_TIMER_INFORMATION
			std::cout << "Elapsed time (ms): " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - startTime).count()) << std::endl;
			std::cout << "Assembling Standard Inverse Mass Two-Normal Operators" << std::endl;
#endif

			addGlobalTwoNormalOperators(res.get());

		}

		res->Finalize();
		return res;
	}

}
