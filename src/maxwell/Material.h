#pragma once

#include <math.h>

namespace maxwell {

class Material {
public:
	Material(const double& ,const  double&);

	const double getPermittivity() const { return epsilon_; }
	const double getPermeability() const { return mu_; }
	const double getImpedance() const { return sqrt(mu_ / epsilon_); }
	const double getConductance() const { return sqrt(epsilon_ / mu_); }

private:
	double epsilon_, mu_;
};

}