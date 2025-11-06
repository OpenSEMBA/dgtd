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

    // Eigenvalues: (2 + 3i), (2 - 3i), (-1)
    // Eigenvectors:
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

    std::complex<double> eigval1(2.0,  3.0);
    std::complex<double> eigval2(2.0, -3.0);
    std::complex<double> eigval3(-1.0, 0.0);

    Eigen::Vector3cd v1_exp = {
        std::complex<double>(1.0,  0.0),
        std::complex<double>(0.0, -1.0),
        std::complex<double>(0.0,  0.0)
    };

    Eigen::Vector3cd v2_exp = {
        std::complex<double>(1.0,  0.0),
        std::complex<double>(0.0,  1.0),
        std::complex<double>(0.0,  0.0)
    };

    Eigen::Vector3cd v3_exp = {
        std::complex<double>(0.0,  0.0),
        std::complex<double>(0.0,  0.0),
        std::complex<double>(1.0,  0.0)
    };

    const double tol = 1e-10;

    //eigval checking
    EXPECT_NEAR(std::real(evals[0]), std::real(eigval1), tol);
    EXPECT_NEAR(std::imag(evals[0]), std::imag(eigval1), tol);

    EXPECT_NEAR(std::real(evals[1]), std::real(eigval2), tol);
    EXPECT_NEAR(std::imag(evals[1]), std::imag(eigval2), tol);

    EXPECT_NEAR(std::real(evals[2]), std::real(eigval3), tol);
    EXPECT_NEAR(std::imag(evals[2]), std::imag(eigval3), tol);

    //eigvec checking
    Eigen::Matrix3cd M;
    M << std::complex<double>(2.0, 0.0), std::complex<double>(-3.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(3.0, 0.0),  std::complex<double>( 2.0, 0.0), std::complex<double>(0.0, 0.0),
         std::complex<double>(0.0, 0.0),  std::complex<double>( 0.0, 0.0), std::complex<double>(-1.0, 0.0);

    auto checkEigenpair = [&](const Eigen::Vector3cd& v, std::complex<double> lambda, const std::string& label)
    {
        Eigen::Vector3cd lhs = M * v;
        Eigen::Vector3cd rhs = lambda * v;
        double err = (lhs - rhs).norm();
        EXPECT_LT(err, tol) << "Eigenvector check failed for " << label << ", residual = " << err;
    };

    checkEigenpair(v1_exp, eigval1, "(2 + 3i)");
    checkEigenpair(v2_exp, eigval2, "(2 - 3i)");
    checkEigenpair(v3_exp, eigval3, "(-1)");
}

TEST_F(SolverExtensionTest, targetIds)
{
    size_t order = 2;
    size_t num_of_segments = 10;
    std::vector<NodeId> ids = buildTargetNodeIds(order, num_of_segments);
    ASSERT_EQ(2, ids[0]);
    ASSERT_EQ(3, ids[1]);
    ASSERT_EQ(33, ids[2]);
    ASSERT_EQ(34, ids[3]);
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

TEST_F(SolverExtensionTest, loadAndUnloadTest)
{

}

}