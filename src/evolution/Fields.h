#pragma once

#include <mfem.hpp>

#include "components/Types.h"

namespace maxwell {

class Fields {
public:
    Fields(mfem::FiniteElementSpace& fes);
    
    mfem::GridFunction& get(const FieldType&, const Direction&);
    const mfem::GridFunction& get(const FieldType&, const Direction&) const;

    mfem::GridFunction& get(const FieldType&);
    
    mfem::Vector& allDOFs() { return allDOFs_; }
    const mfem::Vector& allDOFs() const { return allDOFs_; }

    double getNorml2() const { return allDOFs_.Norml2(); }

private:
    mfem::Vector allDOFs_;
    std::unique_ptr<FiniteElementSpace> global_fes_;
    std::array<mfem::GridFunction, 3> e_, h_;
    mfem::GridFunction e_global_, h_global_;
};
}