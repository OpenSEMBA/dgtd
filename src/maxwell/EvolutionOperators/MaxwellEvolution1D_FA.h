//#pragma once
//
//#include "../mfemExtension/BilinearIntegrators.h"
//
//#include "../Types.h"
//#include "../Model.h"
//#include "../Sources.h"
//#include "../SourcesManager.h"
//#include "../MaxwellEvolutionMethods.h"
//
//#include <Eigen/Dense>
//
//namespace maxwell {
//
//class MaxwellEvolution1D : public mfem::TimeDependentOperator {
//public:
//	static const int numberOfFieldComponents = 2;
//	static const int numberOfMaxDimensions = 1;
//
//	MaxwellEvolution1D(mfem::FiniteElementSpace&, Model&, SourcesManager&, MaxwellEvolOptions&);
//	virtual void Mult(const Vector& x, Vector& y) const;
//	double GetTime() const { return t; }
//	void SetTime(const double time) { t = time; }
//
//	const mfem::FiniteElementSpace& getFES() { return fes_; }
//
//private:
//	std::array<FiniteElementOperator, 2> MS_;
//	std::array<FiniteElementOperator, 2> MF_;
//	std::array<FiniteElementOperator, 2> MP_;
//	std::array<FiniteElementOperator, 2> MT_;
//
//	Eigen::MatrixXd fullMat_;
//
//	mfem::FiniteElementSpace& fes_;
//	Model& model_;
//	SourcesManager& srcmngr_;
//	MaxwellEvolOptions& opts_;
//
//};
//
//}