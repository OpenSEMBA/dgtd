#include <pybind11/pybind11.h>

#include "EvolutionOptions.h"
//#include "SolverOptions.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace maxwell;

PYBIND11_MODULE(solver, m) {
    py::class_<EvolutionOptions>(m, "EvolutionOptions")
        .def(py::init<>())
        .def_readwrite("spectral", &EvolutionOptions::spectral)
        ;
}