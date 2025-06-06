#pragma once

#include <mfem.hpp>

#include "components/Types.h"

namespace maxwell {

class Fields {
public:
    Fields(mfem::ParFiniteElementSpace& fes);
    
    mfem::ParGridFunction& get(const FieldType&, const Direction&);
    const mfem::ParGridFunction& get(const FieldType&, const Direction&) const;

    mfem::ParGridFunction& get(const FieldType&);
    
    mfem::Vector& allDOFs() { return allDOFs_; }
    const mfem::Vector& allDOFs() const { return allDOFs_; }

    double getNorml2() const { return allDOFs_.Norml2(); }

    void updateGlobal();

private:
    mfem::Vector allDOFs_;
    std::unique_ptr<mfem::ParFiniteElementSpace> global_fes_;
    std::array<mfem::ParGridFunction, 3> e_, h_;
    mfem::ParGridFunction e_global_, h_global_;
};
}