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

    H[0] = CampoHx(x, y, z, t);
    H[1] = CampoHy(x, y, z, t);
    H[2] = CampoHz(x, y, z, t);
    E[0] = CampoEx(x, y, z, t);
    E[1] = CampoEy(x, y, z, t);
    E[2] = CampoEz(x, y, z, t);

    // Hopfion::FieldEH CamposEH = <E, H>;
    // Hopfion::FieldEH CamposEH(E, H);
    Hopfion::FieldEH CamposEH = std::make_pair(E, H);

    return CamposEH;

}
