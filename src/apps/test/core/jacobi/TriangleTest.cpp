// OpenSEMBAompu
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.

#include "gtest/gtest.h"
#include "jacobi/Triangle.h"

#include <type_traits>

using namespace SEMBA;
using namespace Math;
using namespace Cudg3d::Jacobi;

class JacobiTriangleTest : public ::testing::Test {
    template<size_t N> friend class Triangle;
protected:
    const Real tol_ = 1e-8;
};

TEST_F(JacobiTriangleTest, LegendreGaussLobatoPoints) {

    {
        static const size_t n = 1;
        Triangle<n> sim;

        std::vector<CVecR2> expected = {
                {-1.000000000000000,  -0.577350269189626},
                { 1.000000000000000,  -0.577350269189626},
                {               0.0,   1.154700538379252} };

        auto computed = sim.getGaussLobattoPoints();

        for (std::size_t i = 0; i < expected.size(); ++i) {
            EXPECT_EQ(expected[i], computed[i]) << "N=" << n;;
        }
    }

//    {
//        static const size_t n = 2;
//        Triangle<n> lin;
//        std::vector<Real> expected = {
//                -1.000000000000000,
//                -0.447213595499958,
//                 0.447213595499958,
//                 1.000000000000000};
//
//        auto computed = lin.getGaussLobattoPoints();
//
//        for (std::size_t i = 0; i < expected.size(); ++i) {
//            EXPECT_FLOAT_EQ(expected[i], computed[i]) << "N=" << n;;
//        }
//    }
//
//    {
//        static const size_t n = 3;
//        Triangle<n> lin;
//        std::vector<Real> expectedLGL = {
//                -1.000000000000000,
//                -0.654653670707977,
//                 0.000000000000000,
//                 0.654653670707977,
//                 1.000000000000000};
//
//        auto computed = lin.getGaussLobattoPoints();
//
//        for (std::size_t i = 0; i < expectedLGL.size(); ++i) {
//            EXPECT_NEAR(expectedLGL[i], computed[i], tol_) << "N=" << n;
//        }
//    }
//
//    {
//        static const size_t n = 4;
//        Triangle<n> lin;
//        std::vector<Real> expected = {
//                -1.000000000000000,
//                -0.765055323929465,
//                -0.285231516480645,
//                 0.285231516480645,
//                 0.765055323929465,
//                 1.000000000000000
//        };
//
//        std::vector<Real> computed = lin.getGaussLobattoPoints();
//
//        for (std::size_t i = 0; i < expected.size(); ++i) {
//            EXPECT_FLOAT_EQ(expected[i], computed[i]) << "N=" << n;
//        }
//    }

}

//TEST_F(JacobiTriangleTest, VandermondeMatrix) {
//    {
//        static const size_t n = 4;
//        Triangle<n> lin;
//        auto r = lin.getGaussLobattoPoints();
//        Matrix::Dynamic<Real> computed = lin.getVandermondeMatrix(r);
//
//        std::vector<Real> expected = {
//            0.707106781186547, -1.224744871391589,  1.581138830084190, -1.870828693386972,  2.121320343559644,
//            0.707106781186547, -0.801783725737273,  0.225876975726313,  0.524890659167824, -0.909137290096990,
//            0.707106781186547, -0.000000000000000, -0.790569415042095,  0.000000000000000,  0.795495128834866,
//            0.707106781186547,  0.801783725737273,  0.225876975726313, -0.524890659167823, -0.909137290096990,
//            0.707106781186547,  1.224744871391589,  1.581138830084190,  1.870828693386972,  2.121320343559644};
//
//        for (size_t i = 0; i < expected.size(); ++i) {
//            EXPECT_NEAR(expected[i], computed.val(i), tol_) << "N=" << n;
//        }
//    }
//}
//
//TEST_F(JacobiTriangleTest, GradVandermondeMatrix) {
//    {
//        static const size_t n = 4;
//        Triangle<n> lin;
//        auto r = lin.getGaussLobattoPoints();
//        Matrix::Dynamic<Real> computed = lin.getGradVandermondeMatrix(r);
//
//        std::vector<Real> expected = {
//                0, 1.224744871391589, -4.743416490252569, 11.224972160321824, -21.213203435596430,
//                0, 1.224744871391589, -3.105295017040594,  3.207134902949095,  -0.000000000000003,
//                0, 1.224744871391589, -0.000000000000000, -2.806243040080456,   0.000000000000001,
//                0, 1.224744871391589,  3.105295017040595,  3.207134902949096,   0.000000000000005,
//                0, 1.224744871391589,  4.743416490252569, 11.224972160321824,  21.213203435596430};
//
//        for (size_t i = 0; i < expected.size(); ++i) {
//            EXPECT_NEAR(expected[i], computed.val(i), tol_) << "N=" << n;
//        }
//    }
//}
//
//TEST_F(JacobiTriangleTest, DifferentiationMatrix) {
//    {
//        static const size_t n = 4;
//        Triangle<n> lin;
//        auto r = lin.getGaussLobattoPoints();
//        Matrix::Dynamic<Real> computed = lin.getDifferentiationMatrix(r);
//
//        std::vector<Real> expected = {
//                -5.000000000000000,  6.756502488724239, -2.666666666666667,  1.410164177942428, -0.500000000000000,
//                -1.240990253030983,  0.000000000000000,  1.745743121887938, -0.763762615825973,  0.259009746969017,
//                 0.375000000000000, -1.336584577695454,  0.000000000000001,  1.336584577695453, -0.375000000000000,
//                -0.259009746969017,  0.763762615825974, -1.745743121887939, -0.000000000000001,  1.240990253030984,
//                 0.500000000000000, -1.410164177942427,  2.666666666666665, -6.756502488724237,  4.999999999999999};
//
//        for (size_t i = 0; i < expected.size(); ++i) {
//            EXPECT_NEAR(expected[i], computed.val(i), tol_) << "N=" << n;
//        }
//    }
//}
//
//TEST_F(JacobiTriangleTest, LIFTMatrix) {
//    {
//        static const size_t n = 3;
//        Triangle<n> lin;
//        Matrix::Dynamic<Real> computed = lin.getLiftMatrix(
//                lin.getGaussLobattoPoints());
//
//        std::vector<Real> expected = {
//                8.000000000000004,  -2.000000000000003,
//               -0.894427190999917,   0.894427190999917,
//                0.894427190999917,  -0.894427190999917,
//               -2.000000000000003,   8.000000000000004};
//
//        EXPECT_EQ(computed.size(), expected.size());
//        for (size_t i = 0; i < computed.size(); ++i) {
//            EXPECT_NEAR(expected[i], computed.val(i), tol_) << "N=" << n;
//        }
//    }
//}
