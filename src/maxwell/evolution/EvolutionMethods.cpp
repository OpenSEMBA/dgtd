#include "EvolutionMethods.h"

namespace maxwell {

using namespace mfem;

using namespace mfemExtension;

InteriorCoefficients intCoeff{
	{FluxType::Centered, { 1.0, 0.0 }},
	{FluxType::Upwind  , { 1.0, 1.0 }}
};

FluxBdrCoefficientsCentered bdrCentCoeff{
	{BdrCond::PEC               ,	{2.0, 0.0}},
	{BdrCond::PMC               ,	{0.0, 2.0}},
	{BdrCond::SMA               ,	{0.0, 0.0}},
	{BdrCond::NONE              ,   {0.0, 0.0}},
};

FluxBdrCoefficientsUpwind bdrUpwindCoeff{
	{BdrCond::PEC               ,	{2.0, 0.0}},
	{BdrCond::PMC               ,	{0.0, 2.0}},
	{BdrCond::SMA               ,	{1.0, 1.0}},
	{BdrCond::NONE              ,   {0.0, 0.0}},
};

FluxSrcCoefficientsCentered srcCentCoeff{
	{BdrCond::TotalFieldInBacked,	{ 1.0, 1.0}},
	{BdrCond::TotalFieldIn      ,	{ 1.0, 1.0}},
	{BdrCond::TotalFieldOut     ,   {-1.0,-1.0}},
};

FluxSrcCoefficientsUpwind srcUpwindCoeff{
	{BdrCond::TotalFieldInBacked,	{ 1.0, 1.0}},
	{BdrCond::TotalFieldIn      ,	{ 1.0, 1.0}},
	{BdrCond::TotalFieldOut     ,   {-1.0,-1.0}},
};

std::map<BdrCond, std::vector<double>> bdrCoeffCheck(const FluxType& ft)
{
	std::map<BdrCond, std::vector<double>> res;
	switch (ft) {
	case(FluxType::Centered):
		res = bdrCentCoeff;
		break;
	case(FluxType::Upwind):
		res = bdrUpwindCoeff;
		break;
	}
	return res;
}

std::map<BdrCond, std::vector<double>> srcCoeffCheck(const FluxType& ft)
{
	std::map<BdrCond, std::vector<double>> res;
	switch (ft) {
	case(FluxType::Centered):
		res = srcCentCoeff;
		break;
	case(FluxType::Upwind):
		res = srcUpwindCoeff;
		break;
	}
	return res;
}

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
	auto res = std::make_unique<BilinearForm>(&fes);
	res->Assemble();
	res->Finalize();
	res->SpMat().Swap(*aux);

	return res;
}

FiniteElementIBFIOperator buildIBFIByMult(
	const BilinearForm& op1,
	const BilinearFormIBFI& op2,
	FiniteElementSpace& fes)
{
	auto aux = mfem::Mult(op1.SpMat(), op2.SpMat());
	auto res = std::make_unique<BilinearFormIBFI>(&fes);
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

	auto MInv = std::make_unique<BilinearForm>(&fes);
	MInv->AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(PWCoeff)));

	MInv->Assemble();
	MInv->Finalize();
	return MInv;
}



FiniteElementOperator buildDerivativeOperator(const Direction& d, FiniteElementSpace& fes)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	if (d >= fes.GetMesh()->Dimension()) {
		res->Assemble();
		res->Finalize();
		return res;
	}

	ConstantCoefficient coeff = (d <= fes.GetMesh()->Dimension()) ? ConstantCoefficient(1.0) : ConstantCoefficient(0.0);
	res->AddDomainIntegrator(
		new DerivativeIntegrator(coeff, d)
	);

	res->Assemble();
	res->Finalize();
	return res;
}



FiniteElementOperator buildFluxOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	res->AddInteriorFaceIntegrator(
		new MaxwellDGTraceJumpIntegrator(dirTerms, intCoeff[opts.fluxType].at(int(FluxType::Centered))));

	for (auto& kv : model.getBoundaryToMarker()) {

		auto c = bdrCoeffCheck(opts.fluxType);
		if (kv.first != BdrCond::SMA) {
			res->AddBdrFaceIntegrator(
				new MaxwellDGTraceJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
		}
		else {
			res->AddBdrFaceIntegrator(
				new MaxwellSMAJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
		}
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	res->AddInteriorFaceIntegrator(
		new MaxwellDGTraceJumpIntegrator(dirTerms, intCoeff[opts.fluxType].at(int(FluxType::Upwind))));

	for (auto& kv : model.getBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		if (kv.first != BdrCond::SMA) {
			res->AddBdrFaceIntegrator(
				new MaxwellDGTraceJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
		}
		else {
			res->AddBdrFaceIntegrator(
				new MaxwellSMAJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
		}
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildZeroNormalOperator(const FieldType& f, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	res->AddInteriorFaceIntegrator(
		new MaxwellDGZeroNormalJumpIntegrator(intCoeff[opts.fluxType].at(int(opts.fluxType))));

	for (auto& kv : model.getBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		res->AddBdrFaceIntegrator(
			new MaxwellDGZeroNormalJumpIntegrator(c[kv.first].at(f)), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildOneNormalOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	res->AddInteriorFaceIntegrator(
		new MaxwellDGOneNormalJumpIntegrator(dirTerms, intCoeff[opts.fluxType].at(int(opts.fluxType))));

	for (auto& kv : model.getBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		res->AddBdrFaceIntegrator(
			new MaxwellDGOneNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildTwoNormalOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);
	res->AddInteriorFaceIntegrator(
		new MaxwellDGTwoNormalJumpIntegrator(dirTerms, intCoeff[opts.fluxType].at(int(opts.fluxType))));

	for (auto& kv : model.getBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		res->AddBdrFaceIntegrator(
			new MaxwellDGTwoNormalJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementOperator buildPenaltyFixOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearForm>(&fes);

	if (model.getInteriorSourceToMarker().size() != 0 && opts.fluxType == FluxType::Upwind) {

		res->AddInteriorFaceIntegrator(
			new MaxwellDGTraceJumpIntegrator(dirTerms, -intCoeff[opts.fluxType].at(int(FluxType::Upwind))));

		for (auto& kv : model.getBoundaryToMarker()) {
			auto c = bdrCoeffCheck(opts.fluxType);
			if (kv.first != BdrCond::SMA) {
				res->AddBdrFaceIntegrator(
					new MaxwellDGTraceJumpIntegrator(dirTerms, -c[kv.first].at(f)), kv.second);
			}
			else {
				res->AddBdrFaceIntegrator(
					new MaxwellSMAJumpIntegrator(dirTerms, -c[kv.first].at(f)), kv.second);
			}
		}
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildFluxIBFIOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		switch (kv.first) {
		case (BdrCond::SMA):
			res->AddInteriorBoundaryFaceIntegrator(
				new MaxwellSMAJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			break;
		default:
			res->AddInteriorBoundaryFaceIntegrator(
				new MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			break;
		}
	}

	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildPenaltyIBFIOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorBoundaryToMarker()) {
		auto c = bdrCoeffCheck(opts.fluxType);
		switch (kv.first) {
		case (BdrCond::SMA):
			res->AddInteriorBoundaryFaceIntegrator(
				new MaxwellSMAJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			break;
		default:
			res->AddInteriorBoundaryFaceIntegrator(
				new MaxwellDGInteriorJumpIntegrator(dirTerms, c[kv.first].at(f)), kv.second);
			break;

		}
	}
	res->Assemble();
	res->Finalize();
	return res;
}

FiniteElementIBFIOperator buildFluxFunctionOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts)
{
	auto res = std::make_unique<BilinearFormIBFI>(&fes);

	for (auto& kv : model.getInteriorSourceToMarker())
	{
		auto c = srcCoeffCheck(opts.fluxType);
		res->AddInteriorBoundaryFaceIntegrator(
			new MaxwellDGFluxTotalFieldIntegrator({ X }, c[kv.first].at(f), 0.5), kv.second
		);
	}

	res->Assemble();
	res->Finalize();
	return res;
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

void calculateEigenvalues(SparseMatrix& mat, Vector& res)
{
	auto denseMat{ mat.ToDenseMatrix() };
	denseMat->Finalize(1);
	denseMat->Eigenvalues(res);
}

double findMaxEigenvalueModulus(const Eigen::VectorXcd& eigvals)
{
	auto res{ 0.0 };
	for (int i = 0; i < eigvals.size(); ++i) {
		auto modulus{ sqrt(pow(eigvals[i].real(),2.0) + pow(eigvals[i].imag(),2.0)) };
		if (modulus <= 1.0 && modulus >= res) {
			res = modulus;
		}
	}
	return res;
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
		throw std::runtime_error("Wrong input in method, check direction or field type vectors.");
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
		throw std::runtime_error("Wrong input in method, check direction or field type vectors.");
	}
	return res;
}

void allocateDenseInEigen1D(DenseMatrix* bilMat, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const double sign)
{
	auto offset = bilMat->Height();
	auto offsetCoeff{ calcOffsetCoeff1D(f) };

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

SparseMatrix toMFEMSparse(const Eigen::SparseMatrix<double>& sp)
{
	SparseMatrix res(int(sp.rows()), int(sp.cols()));
	for (int k = 0; k < sp.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(sp, k); it; ++it)
		{
			res.Set(int(it.row()), int(it.col()), it.value());
		}
	res.Finalize();
	return res;
}

double usePowerMethod(const Eigen::SparseMatrix<double>& global, int iterations)
{
	auto spMat{ toMFEMSparse(global) };
	Vector itVec(int(global.cols()));
	auto mfemSparse{ toMFEMSparse(global) };
	PowerMethod pwrMtd;
	return pwrMtd.EstimateLargestEigenvalue(mfemSparse, itVec, iterations);
}

void performSpectralAnalysis(const FiniteElementSpace& fes, Model& model, const MaxwellEvolOptions& opts)
{
	Array<int> domainAtts(1);
	const auto submAtt{ 501 };
	domainAtts[0] = submAtt;
	auto meshCopy{ fes.GetMesh() };
	auto highestModulus{ 0.0 };
	for (int elem = 0; elem < meshCopy->GetNE(); ++elem) {
		auto preAtt(fes.GetMesh()->GetAttribute(0));
		meshCopy->SetAttribute(elem, domainAtts[0]);
		auto submesh{ SubMesh::CreateFromDomain(*meshCopy,domainAtts) };
		meshCopy->SetAttribute(elem, preAtt);
		auto eigModulus{ findMaxEigenvalueModulus(assembleSubmeshedSpectralOperatorMatrix(submesh, *fes.FEColl(), opts).toDense().eigenvalues()) };
		if (eigModulus >= highestModulus) {
			highestModulus = eigModulus;
			if (highestModulus >= 1.0) {
				std::runtime_error("Modulus of eigenvalue is higher than 1.0 - Unstability.");
			}
		}
	}
}

Eigen::SparseMatrix<double> assembleSubmeshedSpectralOperatorMatrix(Mesh& submesh, const FiniteElementCollection& fec, const MaxwellEvolOptions& opts)
{
	Model submodel(submesh, AttributeToMaterial{}, assignAttToBdrByDimForSpectral(submesh), AttributeToInteriorConditions{});
	FiniteElementSpace subfes(&submesh, &fec);
	Eigen::SparseMatrix<double> local;
	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildDerivativeOperator(y, subfes), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }, -1.0); // MS
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildDerivativeOperator(z, subfes), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildDerivativeOperator(y, subfes), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildDerivativeOperator(z, subfes), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y }, -1.0);

		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildOneNormalOperator(E, { y }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,z }); // MFN
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildOneNormalOperator(E, { z }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,E }, { x,y }, -1.0);
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildOneNormalOperator(H, { y }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,z }, -1.0);
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildOneNormalOperator(H, { z }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,H }, { x,y });

		if (opts.fluxType == FluxType::Upwind) {

			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildZeroNormalOperator(H, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { x }, -1.0); // MP
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildZeroNormalOperator(E, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { x }, -1.0);

			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildTwoNormalOperator(H, { X, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { X,x }); //MPNN
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildTwoNormalOperator(H, { Y, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Y,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, submodel, subfes), *buildTwoNormalOperator(H, { Z, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { H,H }, { Z,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildTwoNormalOperator(E, { X, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { X,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildTwoNormalOperator(E, { Y, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Y,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, submodel, subfes), *buildTwoNormalOperator(E, { Z, x }, submodel, subfes, opts), subfes)->SpMat().ToDenseMatrix(), local, { E,E }, { Z,x });

		}

	}

	return local;
}

AttributeToBoundary assignAttToBdrByDimForSpectral(Mesh& submesh)
{
	switch (submesh.Dimension()) {
	case 1:
		return AttributeToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA} };
	case 2:
		switch (submesh.GetElementType(0)) {
		case Element::TRIANGLE:
			return AttributeToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA } };
		case Element::QUADRILATERAL:
			return AttributeToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
		default:
			std::runtime_error("Incorrect element type for 2D spectral AttToBdr assignation.");
		}
	case 3:
		switch (submesh.GetElementType(0)) {
		case Element::TETRAHEDRON:
			return AttributeToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA} };
		case Element::HEXAHEDRON:
			return AttributeToBoundary{ {1, BdrCond::SMA }, {2, BdrCond::SMA}, {3, BdrCond::SMA }, {4, BdrCond::SMA}, {5, BdrCond::SMA }, {6, BdrCond::SMA} };
		default:
			std::runtime_error("Incorrect element type for 3D spectral AttToBdr assignation.");
		}
	default:
		std::runtime_error("Dimension is incorrect for spectral AttToBdr assignation.");
	}

}
}