#include <vector>
#include <complex>
#include <Eigen/Dense>

#include "Material.h" 
#include "components/Model.h"

namespace maxwell {

struct FitResult {
    Eigen::VectorXcd poles;
    Eigen::VectorXcd residues;
    int N;
    double max_error;
    bool stable;
};

FitResult findOptimalPoles(const SGBCProperties& mat, 
                           double start_freq, 
                           double end_freq, 
                           int num_points = 200, 
                           int max_poles = 10, 
                           double tol = 0.01);

}