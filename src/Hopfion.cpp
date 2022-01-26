#include "Hopfion.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mfem.hpp>

Hopfion::Hopfion(std::size_t pIn, std::size_t qIn)
{
    p = pIn;
    q = qIn;
}

Hopfion::FieldEH Hopfion::evaluate(double time, Vec3 position) const
{
    double x = position[0];
    double y = position[1];
    double z = position[2];

    Hopfion::Vec3 E;
    Hopfion::Vec3 H;

    H[0] = FieldHx(x, y, z, time);
    H[1] = FieldHy(x, y, z, time);
    H[2] = FieldHz(x, y, z, time);
    E[0] = FieldEx(x, y, z, time);
    E[1] = FieldEy(x, y, z, time);
    E[2] = FieldEz(x, y, z, time);

    // Hopfion::FieldEH CamposEH = <E, H>;
    // Hopfion::FieldEH CamposEH(E, H);
    Hopfion::FieldEH CamposEH = std::make_pair(E, H);

    return CamposEH;

}

