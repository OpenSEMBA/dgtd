#include "Spherical.h"

namespace maxwell {

SphericalVector::SphericalVector(const std::vector<double>& p)
{
	assert(p.size() > 1 && p.size() < 4);
	auto t_radius = 0.0;
	for (auto v{ 0 }; v < p.size(); v++) {
		t_radius += std::pow(p[v], 2.0);
	};
	radius = std::sqrt(t_radius);
	p.size() == 2 ? theta = 0.0 : theta = std::acos(p[2] / radius);
	phi = std::atan2(p[1], p[0]);
}

SphericalVector::SphericalVector(const mfem::Vector& p)
{
	assert(p.Size() > 1 && p.Size() < 4);
	auto t_radius = 0.0;
	for (auto v{ 0 }; v < p.Size(); v++) {
		t_radius += std::pow(p[v], 2.0);
	};
	radius = std::sqrt(t_radius);
	p.Size() == 2 ? theta = 0.0 : theta = std::acos(p[2] / radius);
	phi = std::atan2(p[1], p[0]);
}

SphericalVector::SphericalVector(const double r, const double th, const double ph)
{
	radius = r;
	theta = th;
	phi = ph;
}


const std::vector<double> SphericalVector::convertToCartesian() const
{
	std::vector<double> res(3);
	res[0] = radius * std::sin(theta) * std::cos(phi);
	res[1] = radius * std::sin(theta) * std::sin(phi);
	res[2] = radius * std::cos(theta);
	return res;
}

const std::vector<double> SphericalVector::convertSphericalVectorFieldToCartesian(const double Ar, const double At, const double Ap) const
{
	auto cost = std::cos(theta);
	auto sint = std::sin(theta);
	auto cosp = std::cos(phi);
	auto sinp = std::sin(phi);
	std::vector<double> res(3);
	res[0] = cosp * (sint * Ar + cost * At) - sinp * Ap;
	res[1] = sinp * (sint * Ar + cost * At) + cosp * Ap;
	res[2] = cost * Ar - sint * At;
	return res;
}


}