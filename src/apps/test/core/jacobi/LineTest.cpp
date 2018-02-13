// OpenSEMBA
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
#include "jacobi/Line.h"

#include <type_traits>

using namespace SEMBA;
using namespace Math;
using namespace Cudg3d;

template <typename T>
class JacobiLineTypedTest : public ::testing::Test {};

using test_types = ::testing::Types<
    std::integral_constant<std::size_t,2>,
    std::integral_constant<std::size_t,3>,
    std::integral_constant<std::size_t,5>,
    std::integral_constant<std::size_t,8>,
    std::integral_constant<std::size_t,12>>;

TYPED_TEST_CASE(JacobiLineTypedTest, test_types);

TYPED_TEST(JacobiLineTypedTest, BasicOperations) {
    static constexpr std::size_t n = TypeParam::value;
    Jacobi::Line<n> lin;

    Real sum = 0.0;
    std::vector<Real> weights = lin.getWeights();
    for (size_t i = 0; i < weights.size(); ++i) {
        sum += weights[i];
    }
    EXPECT_FLOAT_EQ(1.0, sum);
}

class JacobiLineTest : public ::testing::Test {};

TEST_F(JacobiLineTest, LegendreGaussLobatoPoints) {

    {
        Jacobi::Line<3> lin;
        std::vector<Jacobi::Line<3>::Point> expected = {
                -1.000000000000000,
                -0.447213595499958,
                 0.447213595499958,
                 1.000000000000000};

        auto computed = lin.getPoints();

        for (std::size_t i = 0; i < expected.size(); ++i) {
            EXPECT_FLOAT_EQ(expected[i], computed[i]);
        }
    }

    {
        Jacobi::Line<5> lin;
        std::vector<Jacobi::Line<5>::Point> expected = {
                -1.000000000000000,
                -0.765055323929465,
                -0.285231516480645,
                 0.285231516480645,
                 0.765055323929465,
                 1.000000000000000
        };

        std::vector<Jacobi::Line<5>::Point> computed = lin.getPoints();

        for (std::size_t i = 0; i < expected.size(); ++i) {
            EXPECT_FLOAT_EQ(expected[i], computed[i]);
        }
    }

}
