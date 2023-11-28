#pragma once

#include <mfem.hpp>

#include "components/Types.h"

namespace maxwell {

class Fields {
public:
    Fields(mfem::FiniteElementSpace& fes);

    Fields(const Fields&) = delete;
    Fields(Fields&&) = default;

    Fields& operator=(const Fields&) = delete;
    Fields& operator=(Fields&& gfs) = default;

    ~Fields() = default;
    
    mfem::GridFunction& get(const FieldType&, const Direction&);
    const mfem::GridFunction& get(const FieldType&, const Direction&) const;
    
    mfem::Vector& allDOFs() { return allDOFs_; }
    const mfem::Vector& allDOFs() const { return allDOFs_; }

    double getNorml2() const { return allDOFs_.Norml2(); }

private:
    mfem::Vector allDOFs_;
    std::array<mfem::GridFunction, 3> e_, h_;
};
}