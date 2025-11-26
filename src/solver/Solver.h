#pragma once

#include "ProbesManager.h"
#include "SourcesManager.h"
#include "SolverOptions.h"
#include "SolverExtension.h"

#include "evolution/Fields.h"
#include "evolution/MaxwellEvolution.h"
#include "evolution/GlobalEvolution.h"
#include "evolution/HesthavenEvolution.h"
#include "evolution/HesthavenEvolutionMethods.h"

#include "components/DGOperatorFactory.h"

#include "components/SubMesher.h"
#include "math/PhysicalConstants.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>

namespace maxwell {

std::unique_ptr<ParFiniteElementSpace> buildFiniteElementSpace(ParMesh* m, FiniteElementCollection* fec);
double estimateTimeStep(const Model&, const SolverOptions&, const ParFiniteElementSpace&, const TimeDependentOperator*);
double getMinimumInterNodeDistance(FiniteElementSpace& fes);

class Solver {
public:
    using Vector = mfem::Vector;
    using Position = Vector;
    using ParGridFunction = mfem::ParGridFunction;
    using ODESolver = mfem::ODESolver;
    using ParFields = Fields<ParFiniteElementSpace, ParGridFunction>;
    
    Solver(const Model&, const Probes&, const Sources&, const SolverOptions& = SolverOptions());
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    const ParFields& getFields() const { return fields_; };
    const ParGridFunction& getField(const FieldType& f, const Direction& d) { return fields_.get(f, d); }
    const FieldProbe& getFieldProbe(std::size_t probe) const;
    const PointProbe& getPointProbe(std::size_t probe) const;

    double getTime() const { return time_; }
    double getTimeStep() const { return dt_; }

    Model& getModel() { return model_; }
    ParFiniteElementSpace& getFES() { return *fes_.get(); }

    const mfem::TimeDependentOperator* getFEEvol() const { return evolTDO_.get(); }

    const SolverOptions& getSolverOptions() const { return this->opts_; }

    void run();
    void step();

    void setFinalTime(double final_time) {
        opts_.setFinalTime(final_time);
    }

private:

    SolverOptions opts_;
    Model model_;
    mfem::DG_FECollection fec_;
    std::unique_ptr<mfem::ParFiniteElementSpace> fes_;
    ParFields fields_;
    std::unique_ptr<Device> device_;
    
    SourcesManager sourcesManager_;
    ProbesManager probesManager_;
    
    double time_;
    double dt_;
    std::unique_ptr<ODESolver> odeSolver_;
    
    std::unique_ptr<mfem::TimeDependentOperator> evolTDO_;

    std::unique_ptr<std::vector<SGBCSolver>> sgbcSolvers_;

    void checkOptionsAreValid(const SolverOptions&) const; 
    void assignODESolver();
    std::unique_ptr<TimeDependentOperator> assignEvolutionOperator();
    void initSbcSolvers();
    std::vector<std::pair<NodeId, NodeId>> findSGBCDoFPairs();

    Eigen::SparseMatrix<double> assembleSubmeshedSpectralOperatorMatrix(ParMesh&, const FiniteElementCollection&, const EvolutionOptions&);
    GeomTagToBoundary assignAttToBdrByDimForSpectral(ParMesh&);
    double findMaxEigenvalueModulus(const Eigen::VectorXcd&);
    void performSpectralAnalysis(const ParFiniteElementSpace&, Model&, const EvolutionOptions&);
    void evaluateStabilityByEigenvalueEvolutionFunction(Eigen::VectorXcd& eigenvals, MaxwellEvolution&);
    void writeSimulationStatistics(const Time);
    double calcAverageElementSizeInMesh();
};
}