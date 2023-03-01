#include "MaxwellEvolutionMethods.h"

namespace maxwell {

using namespace mfem;

Eigen::MatrixXd toEigen(const DenseMatrix& mat)
{
	Eigen::MatrixXd res(mat.Width(), mat.Height());
	for (int i = 0; i < mat.Width(); i++) {
		for (int j = 0; j < mat.Height(); j++) {
			res(i, j) = mat.Elem(i, j);
		}
	}
	return res;
}



FiniteElementOperator buildByMult(
	const BilinearForm& op1,
	const BilinearForm& op2,
	FiniteElementSpace& fes)
{
	auto aux = mfem::Mult(op1.SpMat(), op2.SpMat());
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FiniteElementIBFIOperator buildIBFIByMult(
	const BilinearForm& op1,
	const mfemExtension::BilinearFormIBFI& op2,
	FiniteElementSpace& fes)
{
	auto aux = mfem::Mult(op1.SpMat(), op2.SpMat());
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

Vector buildNVector(const Direction& d, const FiniteElementSpace& fes)
{
	const auto dim{ fes.GetMesh()->Dimension() };
	Vector r(dim);
	r = 0.0;
	if (d < dim) {
		r[d] = 1.0;
		return r;
	}
	else {
		return r;
	}
}

FiniteElementOperator buildInverseMassMatrix(const FieldType& f, const Model& model, FiniteElementSpace& fes)
{
	Vector aux{ model.buildPiecewiseArgVector(f) };
	PWConstCoefficient PWCoeff(aux);

	auto MInv = std::make_unique<mfemExtension::BilinearForm>(&fes);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));
		
	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}



FiniteElementOperator buildDerivativeOperator(const Direction& d, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);

	if (d >= fes.GetMesh()->Dimension()) {
		res->Assemble();
		res->Finalize();
		return res;
	}

	ConstantCoefficient coeff(1.0);
	res->AddDomainIntegrator(
		new DerivativeIntegrator(coeff, d)
	);

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);

	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.first);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);

	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient(opts);
		res->AddInteriorFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);

	{
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = interiorPenaltyFluxCoefficient(opts);
		}
		else {
			c = interiorFluxCoefficient();
		}
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		}
		else {
			c = boundaryFluxCoefficient(f, kv.first);
		}
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxJumpOperator(const FieldType& f, const std::vector<Direction>& dirTerms, bool usePenaltyCoefficients, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);

	{
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = interiorPenaltyFluxCoefficient(opts);
		}
		else {
			c = interiorFluxCoefficient();
		}
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker()) {
		FluxCoefficient c;
		if (usePenaltyCoefficients) {
			c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		}
		else {
			c = boundaryFluxCoefficient(f, kv.first);
		}
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildFluxFunctionOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		TFSFOrientationCoefficient c = interiorBoundaryFaceCoefficient(kv.first);
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGFluxTotalFieldIntegrator(dirTerms, c.orient, 1.0), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildPenaltyFunctionOperator(const FieldType& f, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		TFSFOrientationCoefficient c = interiorBoundaryFaceCoefficient(kv.first);
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGPenaltyTotalFieldIntegrator({}, c.orient, 1.0), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildFluxOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);
	auto vec = buildNVector(dirTerms.at(0), fes);

	VectorConstantCoefficient vecCC(vec);
	{
		FluxCoefficient c = interiorFluxCoefficient();
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryFluxCoefficient(f, kv.first);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyOperator1D(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<mfemExtension::BilinearForm>(&fes);
	VectorConstantCoefficient one(Vector({ 1.0 }));
	{
		FluxCoefficient c = interiorPenaltyFluxCoefficient(opts);
		res->AddInteriorFaceIntegrator(new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta));
	}

	for (auto& kv : model.getBoundaryToMarker())
	{
		FluxCoefficient c = boundaryPenaltyFluxCoefficient(f, kv.first, opts);
		res->AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirTerms, c.beta), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildFluxFunctionOperator1D(Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		TFSFOrientationCoefficient c = interiorBoundaryFaceCoefficient(kv.first);
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGFluxTotalFieldIntegrator({ X }, c.orient, 0.5), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildPenaltyFunctionOperator1D(Model& model, FiniteElementSpace& fes)
{
	auto res = std::make_unique<mfemExtension::BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker())
	{
		TFSFOrientationCoefficient c = interiorBoundaryFaceCoefficient(kv.first);
		res->AddInteriorBoundaryFaceIntegrator(
			new mfemExtension::MaxwellDGPenaltyTotalFieldIntegrator({ X }, c.orient, 0.5), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

TFSFOrientationCoefficient interiorBoundaryFaceCoefficient(const BdrCond& bdrCond)
{
	switch (bdrCond)
	{
	case BdrCond::TotalFieldIn:
		return TFSFOrientationCoefficient{  1.0 };
	case BdrCond::TotalFieldOut:
		return TFSFOrientationCoefficient{ -1.0 };
	default:
		throw std::exception("Wrongly defined TotalField attribute");
	}
}

FluxCoefficient interiorFluxCoefficient()
{
	return FluxCoefficient{ 1.0 };
}

FluxCoefficient interiorPenaltyFluxCoefficient(const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0 };
	case FluxType::Upwind:
		return FluxCoefficient{ 1.0 };
	default:
		throw std::exception("No defined FluxType.");
	}
}

FluxCoefficient boundaryFluxCoefficient(const FieldType& f, const BdrCond& bdrC)
{
	switch (bdrC) {
	case BdrCond::PEC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 2.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}
	case BdrCond::PMC:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 2.0 };
		}
	case BdrCond::SMA:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 1.0 };
		case FieldType::H:
			return FluxCoefficient{ 1.0 };
		}
	case BdrCond::TotalFieldIn:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}	
	case BdrCond::NONE:
		switch (f) {
		case FieldType::E:
			return FluxCoefficient{ 0.0 };
		case FieldType::H:
			return FluxCoefficient{ 0.0 };
		}
	default:
		throw std::exception("No defined BdrCond.");
	}
}

FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType& f, const BdrCond& bdrC, const MaxwellEvolOptions& opts)
{
	switch (opts.fluxType) {
	case FluxType::Centered:
		return FluxCoefficient{ 0.0 };
	case FluxType::Upwind:
		switch (bdrC) {
		case BdrCond::PEC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 2.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		case BdrCond::PMC:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 2.0 };
			}
		case BdrCond::SMA:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 1.0 };
			case FieldType::H:
				return FluxCoefficient{ 1.0 };
			}
		case BdrCond::TotalFieldIn:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		case BdrCond::NONE:
			switch (f) {
			case FieldType::E:
				return FluxCoefficient{ 0.0 };
			case FieldType::H:
				return FluxCoefficient{ 0.0 };
			}
		default:
			throw std::exception("No defined BdrCond.");
		}
	default:
		throw std::exception("No defined FluxType.");
	}
}

FieldType altField(const FieldType& f)
{
	switch (f) {
	case E:
		return H;
	case H:
		return E;
	}
	throw std::runtime_error("Invalid field type for altField.");

}

void allocateDenseInEigen(const std::array<FiniteElementOperator, 2>& arr, Eigen::MatrixXd& res, const double sign, bool altField)
{
	int offset = arr[E]->SpMat().ToDenseMatrix()->Height();
	for (int i = 0; i < arr[E]->Height(); ++i) {
		for (int j = 0; j < arr[E]->Width(); ++j) {
			switch (altField) {
			case false:
				res(i, j) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j + offset) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			case true:
				res(i, j + offset) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			}
		}
	}
}

void calculateEigenvalues(const Eigen::MatrixXd& mat, Eigen::VectorXcd& res)
{
	res = mat.eigenvalues();
}

void calculateEigenvalues(const Eigen::SparseMatrix<double>& mat, Eigen::VectorXcd& res)
{
	res = mat.toDense().eigenvalues();
}

void checkEigenvalues(const Eigen::VectorXcd& eigvals)
{
	for (int i = 0; i < eigvals.size(); ++i) {
		//std::cout << eigvals[i].real() << std::endl;
		if (eigvals[i].real() > 1e-2) {
			//throw std::exception("Eigenvalue's real part outside positive tolerance.");
		}
	}
}

void exportSparseToMarketFile(const Eigen::MatrixXd& mat)
{
	Eigen::SparseMatrix<double> sparse = mat.sparseView();
	Eigen::saveMarket(sparse, "SparseMatrix.mtx");
}

Eigen::VectorXd toEigenVector(const Vector& in)
{
	Eigen::VectorXd res;
	res.resize(in.Size());
	for (int i = 0; i < in.Size(); ++i) {
		res(i) = in.Elem(i);
	}
	return res;
}

Vector toMFEMVector(const Eigen::VectorXd& in)
{
	Vector res(int(in.size()));
	for (int i = 0; i < res.Size(); ++i) {
		res(i) = in[i];
	}
	return res;
}

std::vector<int> calcOffsetCoeff1D(const std::vector<FieldType>& f)
{
	std::vector<int> res(2);
	if (f[0] == f[1]) {
		if (f[0] == E) {
			res[0] = 0;
			res[1] = 0;
		}
		else {
			res[0] = 1;
			res[1] = 1;
		}
	}
	else if (f[0] != f[1]) {
		if (f[0] == E) {
			res[0] = 0;
			res[1] = 1;
		}
		else {
			res[0] = 1;
			res[1] = 0;
		}
	}
	else {
		throw std::exception("Wrong input in method, check direction or field type vectors.");
	}
	return res;
}

std::vector<int> calcOffsetCoeff(const std::vector<FieldType>& f, const std::vector<Direction>& d)
{
	std::vector<int> res(2);
	if (d.size() == 1) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = d[0];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = 3 + d[0];
		}
	}
	else if (f[0] == f[1]) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = d[1];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = 3 + d[1];
		}
	}
	else if (f[0] != f[1]) {
		if (f[0] == E) {
			res[0] = d[0];
			res[1] = 3 + d[1];
		}
		else {
			res[0] = 3 + d[0];
			res[1] = d[1];
		}
	}
	else {
		throw std::exception("Wrong input in method, check direction or field type vectors.");
	}
	return res;
}

void allocateDenseInEigen1D(DenseMatrix* bilMat, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const double sign)
{
	auto offset = bilMat->Height();
	auto offsetCoeff{calcOffsetCoeff1D(f)};

	for (int i = 0; i < bilMat->Height(); ++i) {
		for (int j = 0; j < bilMat->Width(); ++j) {
			if (bilMat->Elem(i, j) != 0.0) {
				res.coeffRef(i + offset * offsetCoeff[0], j + offset * offsetCoeff[1]) += sign * bilMat->Elem(i, j);
			}
		}
	}
}

void allocateDenseInEigen(DenseMatrix* bilMat, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign)
{
	auto offset = bilMat->Height();
	auto offsetCoeff{ calcOffsetCoeff(f,d) };

	for (int i = 0; i < bilMat->Height(); ++i) {
		for (int j = 0; j < bilMat->Width(); ++j) {
			if (bilMat->Elem(i, j) != 0.0) {
				res.coeffRef(i + offset * offsetCoeff[0], j + offset * offsetCoeff[1]) += sign * bilMat->Elem(i, j);
			}
		}
	}
}

SparseMatrix toMFEMSparse(Eigen::SparseMatrix<double>& sp)
{
	SparseMatrix res(int(sp.rows()), int(sp.cols()));
	for (int k = 0; k < sp.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(sp, k); it; ++it)
		{
			res.Set(int(it.row()), int(it.col()), it.value());
		}
	return res;
}

double usePowerMethod(Eigen::SparseMatrix<double>& global, int iterations)
{
	auto spMat{ toMFEMSparse(global) };
	Vector itVec(int(global.cols()));
	itVec.Randomize();
	PowerMethod pwrMtd;
	return pwrMtd.EstimateLargestEigenvalue(toMFEMSparse(global), itVec, iterations);
}

}