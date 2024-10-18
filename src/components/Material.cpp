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

Material::Material(double epsilon, double mu, double sigma = 0.0) :
	epsilon_(epsilon),
	mu_(mu),
	sigma_(sigma)
{
	verifyParameters(epsilon, mu, sigma);
}

Material buildVacuumMaterial()
{
	return Material(1.0, 1.0);
}

double Material::getImpedance() const
{
	if (sigma_ == 0.0) {
		return sqrt(mu_ / epsilon_);;
	}
	else {
		throw std::runtime_error("Current implementations does not support impedance calculation for materials with conductivity.");
	}
}

double Material::getAdmitance() const
{
	if (sigma_ == 0.0) {
		return sqrt(epsilon_ / mu_);;
	}
	else {
		throw std::runtime_error("Current implementations does not support admitance calculation for materials with conductivity.");
	}
}

double Material::getSpeedOfWave() const 
{
	if (sigma_ == 0.0) {
		return 1.0 / sqrt(mu_ * epsilon_);
	}
	else {
		throw std::runtime_error("Current implementations does not support wave speed calculation for materials with conductivity.");
	}
}

}