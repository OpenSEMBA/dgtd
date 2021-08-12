#pragma once

#include <array>

class Hopfion {
public:
	typedef std::array<double, 3> Vec3;
	typedef std::pair<Vec3, Vec3> FieldEH;
	
	Hopfion(std::size_t p, std::size_t q);

	FieldEH evaluate(double time, Vec3 position) const;

};