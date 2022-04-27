#pragma once

#include <cmath>

namespace maxwell {

class Material {
public:
	const double epsilon = 1.0;
	const double mu = 1.0;

	Material(const double& epsilon, const double& mu);

	const double getImpedance() const { return sqrt(mu / epsilon); }
	const double getConductance() const { return sqrt(epsilon / mu); }
};

}