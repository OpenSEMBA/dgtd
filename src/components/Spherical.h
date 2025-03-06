#pragma once

#include <vector>
#include <math.h>
#include <cmath>
#include <cassert>
#include <mfem.hpp>

namespace maxwell {

class SphericalVector {
public:

    /*
	* @brief Create a Spherical Vector from three different Cartesian coordinates.
	*/
	SphericalVector(const double x, const double y, const double z);

	/*
	* @brief Create a Spherical Vector from a Cartesian vector.
	* @param[in] p -> Standard Vector of doubles which must have size 2 or 3.
	*/
	SphericalVector(const std::vector<double>& p);

	/*
	* @brief Create a Spherical Vector from a Cartesian MFEM vector.
	* @param[in] p -> MFEM Vector which must have size 2 or 3.
	*/
	SphericalVector(const mfem::Vector& p);

	/*
	* @brief Create a Spherical Vector by inputting its parameters.
	*/
	SphericalVector(const double radius, const double theta, const double phi);

	/*
	* @brief Return a Cartesian coordinate vector of this Spherical vector.
	*/
	const std::vector<double> convertToCartesian() const;

	/*
	* @brief Apply this vector's values to an Spherical Vector Field, then convert it to Cartesian.
	*/
	const std::vector<double> convertSphericalVectorFieldToCartesian(const double Ar, const double At, const double Ap) const;

	double radius;
	double theta;
	double phi;
};

}