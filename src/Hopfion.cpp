#include "Hopfion.h"
#include <iostream>

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

    H[0] = fieldHx(x, y, z, time);
    H[1] = fieldHy(x, y, z, time);
    H[2] = fieldHz(x, y, z, time);
    E[0] = fieldEx(x, y, z, time);
    E[1] = fieldEy(x, y, z, time);
    E[2] = fieldEz(x, y, z, time);

    // Hopfion::FieldEH CamposEH = <E, H>;
    // Hopfion::FieldEH CamposEH(E, H);
    Hopfion::FieldEH CamposEH = std::make_pair(E, H);

    return CamposEH;

}
