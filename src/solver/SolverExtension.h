#include "Solver.h"

namespace maxwell {

struct SGBCTransferMap {

	NodePair src;
	NodePair dst;

};

class SGBCSolver {

	SGBCSolver(const Model&, const Probes&, const Sources&, const SolverOptions&);
	SGBCSolver(const SGBCSolver&) = delete;
	SGBCSolver& operator=(const SGBCSolver&) = delete;
		
private:

    SolverOptions opts_;
    Model model_;
    mfem::DG_FECollection fec_;
    std::unique_ptr<mfem::FiniteElementSpace> fes_;
    Fields fields_;
    std::unique_ptr<Device> device_;

    SourcesManager sourcesManager_;
    ProbesManager probesManager_;

    double time_;
    double dt_;
    std::unique_ptr<ODESolver> odeSolver_{ std::make_unique<mfem::RK4Solver>() };

    std::unique_ptr<mfem::TimeDependentOperator> maxwellEvol_;
    std::vector<std::unique_ptr<mfem::TimeDependentOperator>> sgbcEvol_;

    void initSGBCPreReqs();

    void checkOptionsAreValid(const SolverOptions&) const;
    std::unique_ptr<TimeDependentOperator> assignEvolutionOperator();

    Eigen::SparseMatrix<double> assembleSubmeshedSpectralOperatorMatrix(Mesh&, const FiniteElementCollection&, const EvolutionOptions&);
    GeomTagToBoundary assignAttToBdrByDimForSpectral(Mesh&);
    double findMaxEigenvalueModulus(const Eigen::VectorXcd&);
    void performSpectralAnalysis(const FiniteElementSpace&, Model&, const EvolutionOptions&);
    void evaluateStabilityByEigenvalueEvolutionFunction(Eigen::VectorXcd& eigenvals, MaxwellEvolution&);

};

}