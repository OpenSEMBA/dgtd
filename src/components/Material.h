#pragma once

#include <math.h>
#include <stdexcept>


namespace maxwell {

/**
 * @class Material
 * @brief Electromagnetic material properties in normalized solver units.
 *
 * UNIT CONVENTIONS (CRITICAL FOR SGBC):
 *
 * The Material class stores properties in a normalized system where:
 * - Speed of light c = 1 (normalized)
 * - Vacuum permeability μ₀ = 1 (normalized)
 * - Vacuum permittivity ε₀ = 1 (normalized)
 * - Free space impedance Z₀ = √(μ₀/ε₀) = 1 (normalized)
 *
 * STORED VALUES:
 * - epsilon_: Relative permittivity ε_r (dimensionless, ≥ 1.0)
 *   Physical meaning: ε = ε₀ · ε_r, but ε₀ = 1 in normalized units
 *   Example: For copper at DC, ε_r = 1.0 (to good approximation)
 *
 * - mu_: Relative permeability μ_r (dimensionless, ≥ 1.0)
 *   Physical meaning: μ = μ₀ · μ_r, but μ₀ = 1 in normalized units
 *   Example: For non-magnetic materials, μ_r = 1.0
 *
 * - sigma_: Conductivity normalized as σ_SI × Z₀ [dimensionless in solver]
 *   Conversion from SI units (S/m):
 *     sigma_ = σ_SI [S/m] × 376.73... [Ω] (approximately 377 Ω)
 *   Physical meaning in solver: Loss term in constitutive relation
 *   Example: Copper σ_SI ≈ 5.96e7 S/m → sigma_ ≈ 2.25e10 (normalized)
 *
 * TEMPORAL RELAXATION TIME CALCULATION:
 * For conductive materials, the charge relaxation time is:
 *   τ_SI = ε₀ · ε_r / σ_SI [seconds in SI]
 *
 * To compute in normalized solver units, need to account for the Z₀ normalization:
 *   tau_normalized = epsilon_ / (sigma_ / Z₀_SI) × (ε₀_SI / Z₀_SI)
 *
 * This is CRITICAL for SGBC temporal resolution checking. See:
 *   - driver.cpp::calculateMaximumSourceFrequency() for mesh generation
 *   - SolverExtension.cpp::SGBCWrapper::solve() for temporal criterion
 */
class Material {
public:
	Material(double epsilon, double mu, double sigma_);

	double getPermittivity() const { return epsilon_; }
	double getPermeability() const { return mu_; }
	double getConductivity() const { return sigma_; }
	double getImpedance() const;
	double getAdmittance() const;
	double getSpeedOfWave() const;
private:
	double epsilon_, mu_, sigma_;
};

Material buildVacuumMaterial();

}