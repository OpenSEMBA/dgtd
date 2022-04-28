#include "Material.h"

#include <stdexcept>

namespace maxwell {

Material::Material(double epsilon, double mu) :
	epsilon_(epsilon),
	mu_(mu)
{
	if (epsilon_ < 1.0) {
		throw std::runtime_error("Invalid epsilon.");
	}
	if (mu_ < 1.0) {
		throw std::runtime_error("Invalid mu.");
	}
}

}