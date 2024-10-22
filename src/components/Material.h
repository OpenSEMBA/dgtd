#pragma once

#include <math.h>

namespace maxwell {

class Material {
public:
	Material(double epsilon, double mu, double sigma_);

	double getPermittivity() const { return epsilon_; }
	double getPermeability() const { return mu_; }
	double getConductivity() const { return sigma_; }
	double getImpedance() const;
	double getAdmitance() const;
	double getSpeedOfWave() const;
private:
	double epsilon_, mu_, sigma_;
};

Material buildVacuumMaterial();

}