#include <iostream>
#include <cmath>
#include <iomanip>

#include <Eigen/Eigenvalues>

#include "SIBCPoleResidueFitting.h"
#include "math/PhysicalConstants.h"

using Complex = std::complex<double>;
using VectorCx = Eigen::VectorXcd;
using VectorCd = Eigen::VectorXd;
using MatrixCx = Eigen::MatrixXcd;

namespace maxwell{


struct NormalizedFreqs {
    VectorCx s;
    VectorCx s_norm;
    double scale;
};

struct LinearSystem {
    MatrixCx A;
    VectorCx b;
};

static VectorCx generatePhysicsTarget(const SGBCProperties& props, 
                                      const std::vector<double>& freqs_vec) 
{
    auto num_pts = freqs_vec.size();
    VectorCx Z_target(num_pts);
    
    const double eps0 = physicalConstants::vacuumPermittivity_SI;
    const double mu0  = physicalConstants::vacuumPermeability_SI;

    double eps_r = props.layers.front().material.getPermittivity();
    double mu_r  = props.layers.front().material.getPermeability();
    double sigma = props.layers.front().material.getConductivity();
    double thickness = props.totalWidth();

    for(int i = 0; i < num_pts; ++i) {
        double w = 2.0 * M_PI * freqs_vec[i];
        
        Complex eps_c = (eps0 * eps_r) - Complex(0.0, 1.0) * (sigma / w);
        Complex mu_c  = mu0 * mu_r;

        Complex eta = std::sqrt(mu_c / eps_c);
        Complex gamma = Complex(0.0, 1.0) * w * std::sqrt(mu_c * eps_c);

        Z_target[i] = eta / std::tanh(gamma * thickness);
    }
    
    return Z_target;
}

static NormalizedFreqs scaleFrequencies(const std::vector<double>& freqs_vec) 
{
    long num_pts = freqs_vec.size();
    VectorCx s(num_pts);
    double max_mag = 0.0;

    for(int i = 0; i < num_pts; ++i) {
        s[i] = Complex(0.0, 1.0) * (2.0 * M_PI * freqs_vec[i]);
        if(std::abs(s[i]) > max_mag) max_mag = std::abs(s[i]);
    }

    return {s, s / max_mag, max_mag};
}

static VectorCd calculate_weights(const VectorCx& Z_target) 
{
    long num_pts = Z_target.size();
    bool relative_mode = (std::abs(Z_target[0]) < std::abs(Z_target[num_pts - 1]));
    
    VectorCd weights(num_pts);
    for(int i = 0; i < num_pts; ++i) {
        if(relative_mode) {
            weights[i] = 1.0 / (std::abs(Z_target[i]) + 1e-15);
        } else {
            weights[i] = 1.0;
        }
    }
    return weights;
}

static LinearSystem build_system(const NormalizedFreqs& freqs, 
                                 const VectorCx& Z_target, 
                                 const VectorCd& weights, 
                                 int order) 
{
    long num_pts = freqs.s.size();
    long num_vars = 2 * order + 1;
    MatrixCx A(num_pts, num_vars);
    VectorCx b(num_pts);

    for(int k = 0; k < num_pts; ++k) {
        std::complex<double> w = weights[k];
        std::complex<double> z = Z_target[k];
        std::complex<double> s_val = freqs.s_norm[k];

        // Fill B(s) terms: 1, s, s^2 ...
        std::complex<double> p = 1.0;
        for(int j = 0; j <= order; ++j) {
            A(k, j) = w * p; 
            p *= s_val;
        }
        
        p = s_val; 
        for(int j = 0; j < order; ++j) {
            A(k, order + 1 + j) = -w * z * p;
            p *= s_val;
        }

        b(k) = w * z;
    }

    return {A, b};
}

static VectorCx extractRoots(const VectorCx& coeffs_descending) 
{
    long n = coeffs_descending.size(); 
    if (n < 1) return VectorCx(0);

    MatrixCx C = MatrixCx::Zero(n, n);
    for (long i = 1; i < n; ++i) C(i, i - 1) = 1.0;
    
    for (long i = 0; i < n; ++i) C(i, n - 1) = -coeffs_descending(n - 1 - i);

    Eigen::ComplexEigenSolver<MatrixCx> solver(C, false);
    return solver.eigenvalues();
}

static VectorCx solveForPoles(const LinearSystem& sys, double scale, int order) {

    VectorCx x = sys.A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(sys.b);
    VectorCx a_poly(order);

    VectorCx a_full(order + 1);
    a_full(0) = 1.0;
    
    for(int i = 0; i < order; ++i) {
        a_full(i + 1) = x(order + 1 + i);
    }

    for(int i = 0; i <= order; ++i) {
        a_full(i) /= std::pow(scale, i);
    }

    std::complex<double> a_N = a_full(order);
    a_full /= a_N;

    VectorCx coeffs_for_roots(order);
    for(int i = 0; i < order; ++i) {
        coeffs_for_roots(i) = a_full(order - 1 - i);
    }
    
    return extractRoots(coeffs_for_roots);
}

static VectorCx calcResidues(const NormalizedFreqs& freqs, 
                                           const VectorCx& Z_target, 
                                           const VectorCx& poles) {
    long num_pts = freqs.s.size();
    long num_poles = poles.size();

    MatrixCx R_mat(num_pts, num_poles);
    for(int k = 0; k < num_pts; ++k) {
        for(int j = 0; j < num_poles; ++j) {
            R_mat(k, j) = 1.0 / (freqs.s[k] - poles[j]);
        }
    }

    return R_mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Z_target);
}

static FitResult checkQuality(const NormalizedFreqs& freqs,
                               const VectorCx& Z_target,
                               const VectorCx& poles,
                               const VectorCx& residues,
                               int order) 
{
    
    bool stable = true;
    for(int i = 0; i < poles.size(); ++i) {
        if(poles[i].real() > 1e-6) stable = false;
    }

    double max_err = 0.0;
    long num_pts = freqs.s.size();
    
    for(int k = 0; k < num_pts; ++k) {
        std::complex<double> z_fit = 0.0;
        for(int j = 0; j < order; ++j) {
            z_fit += residues[j] / (freqs.s[k] - poles[j]);
        }
        
        double err = std::abs(z_fit - Z_target[k]) / (std::abs(Z_target[k]) + 1e-20);
        if(err > max_err) max_err = err;
    }

    return {poles, residues, order, max_err, stable};
}

FitResult calcPolesAndResidues(const std::vector<double>& freqs_vec, const VectorCx& Z_target, int order)
{
    
    auto freqs = scaleFrequencies(freqs_vec);
    auto weights = calculate_weights(Z_target);

    auto system = build_system(freqs, Z_target, weights, order);

    auto poles = solveForPoles(system, freqs.scale, order);
    auto residues = calcResidues(freqs, Z_target, poles);

    return checkQuality(freqs, Z_target, poles, residues, order);
}

FitResult findOptimalPoles(const SGBCProperties& props, 
                           double start_freq, 
                           double end_freq, 
                           int num_points, 
                           int max_poles, 
                           double tol) 
{

    std::vector<double> freqs(num_points);
    double step = std::pow(end_freq / start_freq, 1.0 / (num_points - 1));
    freqs[0] = start_freq;
    for(int i = 1; i < num_points; ++i) freqs[i] = freqs[i-1] * step;

    VectorCx Z_target = generatePhysicsTarget(props, freqs);

    FitResult res;
    res.max_error = 1e9;
    res.N = 0;
    res.stable = false;

    for(int N = 1; N <= max_poles; ++N) {
        
        FitResult temp_res = calcPolesAndResidues(freqs, Z_target, N);

        if(!temp_res.stable){
            continue;
        }
            
        if(temp_res.max_error < res.max_error) {
            res = temp_res;
        }

        if(temp_res.max_error < tol) {
            return temp_res; 
        }
    }

    return res;
}

}