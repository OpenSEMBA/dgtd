#include <pybind11/pybind11.h>

#include "Solver.h"

namespace py = pybind11;

using namespace maxwell;

PYBIND11_MODULE(maxwell_solver, m) {

    py::enum_<FluxType>(m, "FluxType")
        .value("Centered", FluxType::Centered)
        .value("Upwind", FluxType::Upwind)
        ;

    py::class_<EvolutionOptions>(m, "EvolutionOptions")
        .def(py::init<>())
        .def_readwrite("fluxType",    &EvolutionOptions::fluxType)
        .def_readwrite("spectral",    &EvolutionOptions::spectral)
        .def_readwrite("eigenvals",   &EvolutionOptions::eigenvals)
        .def_readwrite("marketFile",  &EvolutionOptions::marketFile)
        .def_readwrite("powerMethod", &EvolutionOptions::powerMethod)
        ;

    py::class_<SolverOptions>(m, "SolverOptions")
        .def(py::init<>())
        .def_readwrite("order",      &SolverOptions::order)
        .def_readwrite("dt",         &SolverOptions::dt)
        .def_readwrite("t_final",    &SolverOptions::t_final)
        .def_readwrite("CFL",        &SolverOptions::CFL)
        .def_readwrite("evolutionOperatorOptions", 
                                     &SolverOptions::evolutionOperatorOptions)
        ;

    //py::class_<ProblemDescription>(m, "ProblemDescription")
    //    .def(py::init<const Model&, const Probes&, const Sources&>())
    //    ;

    //py::class_<maxwell::Solver>(m, "Solver")
    //    .def(py::init<const ProblemDescription&, const SolverOptions&>())
    //    ;
}