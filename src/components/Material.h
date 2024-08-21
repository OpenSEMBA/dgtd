#pragma once

#include <math.h>

namespace maxwell {

class Material {
public:
	Material(double epsilon, double mu);
	Material(double epsilon, double mu, double sigma_);

	double getPermittivity() const { return epsilon_; }
	double getPermeability() const { return mu_; }
	double getConductivity() const { return sigma_; }
	double getImpedance() const { return sqrt(mu_ / epsilon_); }
	double getAdmitance() const { return sqrt(epsilon_ / mu_); }
	double getSpeedOfLight() const { return 1 / sqrt(mu_ * epsilon_); }
private:
	double epsilon_, mu_, sigma_;
};

Material buildVacuumMaterial();

}