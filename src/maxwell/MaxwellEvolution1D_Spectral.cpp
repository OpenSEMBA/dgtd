#include "MaxwellEvolution1D_Spectral.h"
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

void allocateDenseInEigen1D(const std::array<FiniteElementOperator, 2>& arr, Eigen::MatrixXd& res, const double sign = 1.0, bool altField = false)
{
	int offset = arr[E]->SpMat().ToDenseMatrix()->Height();
	for (int i = 0; i < arr[E]->Height(); ++i) {
		for (int j = 0; j < arr[E]->Width(); ++j) {
			switch (altField) {
			case false:
				res(i         , j         ) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j + offset) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			case true:
				res(i         , j + offset) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j         ) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			}
		}
	}
}

void allocateDenseInEigen1D(const std::array<FiniteElementIBFIOperator, 2>& arr, Eigen::MatrixXd& res, const double sign = 1.0, bool altField = false)
{
	int offset = arr[E]->SpMat().ToDenseMatrix()->Height();
	for (int i = 0; i < arr[E]->Height(); ++i) {
		for (int j = 0; j < arr[E]->Width(); ++j) {
			switch (altField) {
			case false:
				res(i         , j         ) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j + offset) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			case true:
				res(i         , j + offset) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j         ) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			}
		}
	}
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

void calculateEigenvalues(const Eigen::MatrixXd& mat, Eigen::VectorXcd& res)
{
	res = mat.eigenvalues();
}

void checkEigenvalues(const Eigen::VectorXcd eigvals)
{
	for (int i = 0; i < eigvals.size(); ++i) {
		if (eigvals[i].real() > 1e-10) {
			throw std::exception("Eigenvalue's real part outside positive tolerance.");
		}
	}
}

void exportSparseToMarketFile(const Eigen::MatrixXd& mat)
{
	Eigen::SparseMatrix<double> sparse = mat.sparseView();
	Eigen::saveMarket(sparse, "SparseMatrix.mtx");
}

MaxwellEvolution1D_Spectral::MaxwellEvolution1D_Spectral(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
	std::array<FiniteElementOperator, 2> MS, MF, MP;
	std::array<FiniteElementIBFIOperator, 2> MBF, MBP;
	for (auto f : {E, H}) {
		const auto f2{ altField(f) };
		MS[f] = buildByMult		 (*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(X, fes_), fes_);
		MF[f] = buildByMult		 (*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator1D(f2, {X}, model_, fes_), fes_);
		MBF[f] = buildIBFIByMult	 (*buildInverseMassMatrix(f, model_, fes_), *buildFluxFunctionOperator1D(model_, fes_), fes_);
		if (opts_.fluxType == FluxType::Upwind) {
			MP[f] = buildByMult	 (*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyOperator1D(f, {}, model_, fes_, opts_), fes_);
			MBP[f] = buildIBFIByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyFunctionOperator1D(model_, fes_), fes_);
		}
	}

	global_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
				   numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	global_.setZero();
	allocateDenseInEigen1D(MS, global_, -1.0, true);
	allocateDenseInEigen1D(MF, global_,  1.0, true);

	if (opts_.fluxType == FluxType::Upwind) {
		allocateDenseInEigen1D(MP, global_, -1.0);
	}

	forcing_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	forcing_.setZero();

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<PlaneWave*>(source.get())) {
			allocateDenseInEigen1D(MBF, forcing_);
			if (opts_.fluxType == FluxType::Upwind) {
				allocateDenseInEigen1D(MBP, forcing_);
			}
		}
	}

	calculateEigenvalues(global_, eigenvals_);
	checkEigenvalues(eigenvals_);
	exportSparseToMarketFile(global_);

}


void MaxwellEvolution1D_Spectral::Mult(const Vector& in, Vector& out) const
{
	Eigen::VectorXd fieldsOld{ toEigenVector(in) };
	Eigen::VectorXd fieldsNew{ global_ * fieldsOld };

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<PlaneWave*>(source.get())) {
			GridFunction eFunc(srcmngr_.evalTotalField(GetTime()));
			GridFunction hFunc(srcmngr_.evalTotalField(GetTime()));

			Eigen::VectorXd forcVecsOld;
			forcVecsOld.resize(2 * eFunc.Size());
			for (int i = 0; i < eFunc.Size(); ++i) {
				forcVecsOld(i)                = eFunc.Elem(i);
				forcVecsOld(i + eFunc.Size()) = hFunc.Elem(i);
			}
			fieldsNew += forcing_ * forcVecsOld;
		}
	}
	out = toMFEMVector(fieldsNew);
}

}

