#include "SquareResonantCavity2D.h"
#include "mfem.hpp"
#include "driver.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace mfem;

SquareResonantCavity2D::SquareResonantCavity2D()
{
	// Pensar esto mejor: Como trabajo con el mallado dado por mfem
    // --> Una opcion es obtener y calcular los valores del mallado desde aquí
    // --> Otra puede ser que Driver::Run devuelva todos los valores que necesito como matriz 
    //      similar a como se procedió con el Test del hopfion.
    // 
    // 
	// tiempo a usar = time;
	// matriz de posiciones de mfem = position;

    // Si se toma la primera opción sería Algo así?

    /*
    * 
    std::cout.precision(Driver.Options.precision);

    Device device(Driver.Options.device_config);

    Mesh mesh(Driver.Options.mesh_file, 1, 1);

    int dim = mesh.Dimension();

    for (int lev = 0; lev < Driver.Options.ref_levels; lev++)
    {
        mesh.UniformRefinement();
    }

    // De aqui obtener las posiciones del grid (X,Y,Z)
    mesh.GetBoundingBox(meshBoundingBoxMin, meshBoundingBoxMax, std::max(Driver.Options.order, 1));
    
    // De aquí obtener el tiempo
    t = Driver.Options.t_final 
    
    */
}

SquareResonantCavity2D::FieldEH SquareResonantCavity2D::TheoField(double time, Vec3 position) const
{
	// Calculo teorico del campo EM en la cavidad resonate a un tiempo "time" y en la posicion "position"

    double x = position[0];
    double y = position[1];
    double z = position[2];

    // Expresion de los campos

    Hopfion::Vec3 E;
    Hopfion::Vec3 H;

    H[0] = theofieldHx(x, y, z, time);
    H[1] = theofieldHy(x, y, z, time);
    H[2] = theofieldHz(x, y, z, time);
    E[0] = theofieldEx(x, y, z, time);
    E[1] = theofieldEy(x, y, z, time);
    E[2] = theofieldEz(x, y, z, time);

    // Hopfion::FieldEH CamposEH = <E, H>;
    // Hopfion::FieldEH CamposEH(E, H);
    SquareResonantCavity2D::FieldEH CamposEH = std::make_pair(E, H);

    return CamposEH;

}

SquareResonantCavity2D::FieldEH SquareResonantCavity2D::CompuField(double time) const
{
	// Calculo numerico del campo en la cavidad resonate hasta un tiempo t
	// Puedo intentar que esta devuelva el tiempo, el mallado y los campos.
	//   con el valor del tiempo y del mallado genero los valores teorico.
    // Pero entonces no puede ser FieldEH, hay que crear una clase nueva
    // 
    // Driver.Options.paraview = false tiene que serlo para devolver los datos como quiero
    // E, H = Driver.run() o así del estilo.
    // return t, x, y, z, E[0], E[1], E[2], H[0], H[1], H[2].
}