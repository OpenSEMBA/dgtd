#include <gtest/gtest.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "solver/SolverExtension.h"
#include "solver/Solver.h"
#include "cases/CasesFunctions.h"
#include "driver/driver.h"
#include "TestUtils.h"


namespace maxwell{

using namespace mfem;
using json = nlohmann::json;


class SolverExtensionTest : public ::testing::Test {
protected:

    static json parseJSONfile(const std::string& json_file)
    {
	    std::ifstream test_file(json_file);
	    return json::parse(test_file);
    }

    static bool checkIfInsideTargetIds(const Nodes& target_ids, const size_t localdofs, const size_t totaldofs, const size_t c){
        bool check = false;
        for (int id : target_ids)
        {
            for (int offset = 0; offset < totaldofs; offset += localdofs)
            {
                if (c == id + offset)
                {
                    check = true;
                    break;
                }
            }
            if (check){
                break;
            }
        }
        return check;
    }

};

TEST_F(SolverExtensionTest, isCorrect_SGBC_Properties)
{
    // auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
    // auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SGBC_Test"), true);
    // for (auto props : global_solver.getSolverOptions().sbc_props){
    //     if (props.phys_tag == 15){
    //         ASSERT_EQ(1e-5, props.material_width);
    //         ASSERT_EQ(15, props.num_of_segments);
    //         ASSERT_EQ(5, props.order);
    //         ASSERT_EQ(1e5, props.material.getConductivity());
    //     }
    //     else if (props.phys_tag == 18){
    //         ASSERT_EQ(1e-8, props.material_width);
    //         ASSERT_EQ(18, props.num_of_segments);
    //         ASSERT_EQ(8, props.order);
    //         ASSERT_EQ(1e8, props.material.getConductivity()); // High only so it ends in 8, for 'same-number' easy checks.
    //     }
    // }
    ASSERT_TRUE(false); // While reworking phys_tag
}

TEST_F(SolverExtensionTest, checkEigenSolverSolutions)
{
    // A = [ 2  -3   0 ]
    //     [ 3   2   0 ]
    //     [ 0   0  -1 ]
    //
    // Eigenvalues: (2 + 3i), (2 - 3i), (-1)
    // Eigenvectors (one choice):
    //   For (2 + 3i): [1, -i, 0]^T
    //   For (2 - 3i): [1,  i, 0]^T
    //   For (-1):     [0,  0, 1]^T

    SparseMatrix A(3, 3);
    A.Add(0, 0,  2.0);
    A.Add(0, 1, -3.0);
    A.Add(1, 0,  3.0);
    A.Add(1, 1,  2.0);
    A.Add(2, 2, -1.0);
    A.Finalize();

    auto es = applyEigenSolverOnGlobalOperator(A);

    Eigen::VectorXcd evals = es.eigenvalues();
    Eigen::MatrixXcd evecs = es.eigenvectors();

    ASSERT_EQ(3, evals.size());
    ASSERT_EQ(3, evecs.cols());
    ASSERT_EQ(3, evecs.rows());

    const double tol = 1e-10;

    // Expected eigenvalues
    std::vector<std::complex<double>> expected_evals = {
        std::complex<double>(2.0,  3.0),
        std::complex<double>(2.0, -3.0),
        std::complex<double>(-1.0, 0.0)
    };

    Eigen::Vector3cd v1_exp; // (2 + 3i)
    v1_exp << std::complex<double>(1.0, 0.0),
              std::complex<double>(0.0,-1.0),
              std::complex<double>(0.0, 0.0);

    Eigen::Vector3cd v2_exp; // (2 - 3i)
    v2_exp << std::complex<double>(1.0, 0.0),
              std::complex<double>(0.0, 1.0),
              std::complex<double>(0.0, 0.0);

    Eigen::Vector3cd v3_exp; // (-1)
    v3_exp << std::complex<double>(0.0, 0.0),
              std::complex<double>(0.0, 0.0),
              std::complex<double>(1.0, 0.0);

    std::vector<Eigen::Vector3cd> expected_evecs = { v1_exp, v2_exp, v3_exp };

    auto findBestMatchIndex = [&](const std::complex<double>& target, const Eigen::VectorXcd& vals, const std::vector<int>& used) -> int
    {
        double bestDist = std::numeric_limits<double>::infinity();
        int bestIdx = -1;
        for (int i = 0; i < vals.size(); ++i)
        {
            if (std::find(used.begin(), used.end(), i) != used.end()) continue; // skip already used
            double dist = std::abs(vals[i] - target);
            if (dist < bestDist)
            {
                bestDist = dist;
                bestIdx = i;
            }
        }
        return bestIdx;
    };

    auto compareEigenvectors = [&](const Eigen::Vector3cd& v_computed,
                                   const Eigen::Vector3cd& v_expected,
                                   const std::string& label)
    {
        Eigen::Vector3cd vc = v_computed.normalized();
        Eigen::Vector3cd ve = v_expected.normalized();

        std::complex<double> phase = (ve.adjoint() * vc)(0,0);
        if (std::abs(phase) > 0.0)
        {
            vc /= (phase / std::abs(phase));
        }

        double diff = (vc - ve).norm();
        EXPECT_LT(diff, tol) << "Eigenvector mismatch for " << label << " (difference = " << diff << ")";
    };

    std::vector<int> used_indices;

    for (size_t k = 0; k < expected_evals.size(); ++k)
    {
        const auto& exp_ev = expected_evals[k];
        int idx = findBestMatchIndex(exp_ev, evals, used_indices);
        ASSERT_NE(idx, -1) << "Could not find a matching computed eigenvalue for expected eigenvalue index " << k;

        used_indices.push_back(idx);

        EXPECT_NEAR(std::real(evals[idx]), std::real(exp_ev), tol) << "Eigenvalue real part mismatch for expected index " << k;
        EXPECT_NEAR(std::imag(evals[idx]), std::imag(exp_ev), tol) << "Eigenvalue imag part mismatch for expected index " << k;

        compareEigenvectors(evecs.col(idx), expected_evecs[k], std::to_string(k) + " (expected " + std::to_string(std::real(exp_ev)) + (std::imag(exp_ev) >= 0 ? "+" : "") + std::to_string(std::imag(exp_ev)) + "i)");
    }

    // Additionally, perform a sanity check that A * v = lambda * v holds for the matched pairs
    // (This is redundant if eigenvalues/eigenvectors matched above, but kept as an extra safeguard.)
    Eigen::Matrix3cd M;
    M << std::complex<double>(2.0, 0.0), std::complex<double>(-3.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(3.0, 0.0),  std::complex<double>( 2.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0),  std::complex<double>( 0.0, 0.0), std::complex<double>(-1.0, 0.0);

    for (size_t k = 0; k < expected_evals.size(); ++k)
    {
        int matched_idx = used_indices[k];
        Eigen::Vector3cd lhs = M * evecs.col(matched_idx);
        Eigen::Vector3cd rhs = evals[matched_idx] * evecs.col(matched_idx);
        double err = (lhs - rhs).norm();
        EXPECT_LT(err, 1e-8) << "Sanity Av - lambda v residual too large for matched index " << matched_idx << " (residual = " << err << ")";
    }
}

TEST_F(SolverExtensionTest, targetIds)
{
    size_t order = 1;
    size_t num_of_segments = 10;
    std::vector<NodeId> ids = buildTargetNodeIds(order, num_of_segments);
    ASSERT_EQ(1, ids[0]);
    ASSERT_EQ(2, ids[1]);
    ASSERT_EQ(21, ids[2]);
    ASSERT_EQ(22, ids[3]);
}

TEST_F(SolverExtensionTest, buildTest)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SGBC_Test"), true);
    std::pair<GlobalId,GlobalId> pair({2,3});
    for (const auto prop : global_solver.getSolverOptions().sbc_props){
        ASSERT_NO_THROW(SGBCSolver(&prop, pair));
    }
}

TEST_F(SolverExtensionTest, loadAndUnloadSGBCValuesTest)
{
    Material mat(1.0, 1.0, 1e4);
    SGBCProperties props(mat);

    std::pair<GlobalId, GlobalId> ids(1, 2);

    SGBCSolver solver(&props, ids);
    SGBCNodalFields nodal;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            nodal.at(f).at(d).first = double(f + d + 1);
            nodal.at(f).at(d).second = 2.0 * nodal.at(f).at(d).first;
        }
    }

    solver.setSGBCFieldValues(nodal);
    auto returned = solver.getSGBCFieldValues();

    double tol = 1e-3;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            EXPECT_NEAR(nodal.at(f).at(d).first, returned.at(f).at(d).first, tol);
            EXPECT_NEAR(nodal.at(f).at(d).second, returned.at(f).at(d).second, tol);
        }
    }

}

TEST_F(SolverExtensionTest, loadAndUnloadFullStateTest)
{
    Material mat(1.0, 1.0, 1e4);
    SGBCProperties props(mat);
    std::pair<GlobalId, GlobalId> ids(1, 2);
    int ndofs = (props.num_of_segments + 2) * (props.order + 1);

    SGBCSolver solver(&props, ids);
    FullNodalFields nodal;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            nodal.at(f).at(d).SetSize(ndofs);
            for (auto v = 0; v < ndofs; v++){
                nodal.at(f).at(d)[v] = double(f + d + v + 1);
            }
        }
    }

    solver.setFullNodalState(nodal);
    auto returned = solver.getFullNodalState();

    double tol = 1e-3;
    for (auto f : {E, H}){
        for (auto d : {X, Y, Z}){
            for (auto v = 0; v < ndofs; v++){
                EXPECT_NEAR(nodal.at(f).at(d)[v], returned.at(f).at(d)[v], tol);
                EXPECT_NEAR(nodal.at(f).at(d)[v], returned.at(f).at(d)[v], tol);
            }        
        }
    }
}

}