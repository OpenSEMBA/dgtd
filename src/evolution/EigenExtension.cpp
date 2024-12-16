#include "evolution/HesthavenEvolutionTools.h"

namespace maxwell {

    struct MatrixCompare {

        bool operator()(const Eigen::MatrixXd& lhs, const Eigen::MatrixXd& rhs) const {
            
            if (lhs.norm() != rhs.norm()) {
                return lhs.norm() < rhs.norm();  
            }
            
            for (int i{ 0 }; i < lhs.rows(); i++) {
                for (int j{ 0 }; j < lhs.cols(); j++) {
                    if (lhs(i, j) != rhs(i, j)) {
                        return lhs(i, j) < rhs(i, j);  
                    }
                }
            }
            return false;
        }

    };

}