#include "Material.h"

#include <stdexcept>

namespace maxwell {

void verifyParameters(double epsilon, double mu, double sigma) 
{
	if (epsilon < 1.0) {
		throw std::runtime_error("Permittivity under 1.0 not allowed.");
	}
	if (mu < 1.0) {
		throw std::runtime_error("Permeability under 1.0 not allowed.");
	}
	if (sigma < 0.0) {
		throw std::runtime_error("Conductivity under 0.0 not allowed.");
	}
}

Material::Material(double epsilon, double mu, double sigma) :
	epsilon_(epsilon),
	mu_(mu),
	sigma_(sigma)
{
	verifyParameters(epsilon, mu, sigma);
}

Material buildVacuumMaterial()
{
	return Material(1.0, 1.0, 0.0);
}

}