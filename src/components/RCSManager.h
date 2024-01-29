#pragma once

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>

#include <components/Probes.h>


namespace maxwell {

using namespace mfem;
using Rho = double;
using Phi = double;
using SphericalAngles = std::vector<std::pair<Rho, Phi>>;
using Nedelec_XY = GridFunction;
using H1_Z = GridFunction;

class RCSManager {
public:

	RCSManager(const std::string& path, const std::vector<double>& frenquency, const SphericalAngles& angle);

private:

	void performRCS2DCalculations(GridFunction& Ax, GridFunction& Ay, GridFunction& Bz, const double frequency, const std::pair<Rho,Phi>&);

	Mesh m_;
	std::string basePath_;

};

}
