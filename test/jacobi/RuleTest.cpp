#include "gtest/gtest.h"
#include "jacobi/Rule.h"

#include <type_traits>
#include <array>

using namespace SEMBA;
using namespace Math;
using namespace dgtd::jacobi;

class RuleTest : public ::testing::Test {};

TEST_F(RuleTest, QuadraturePointsAndWeights) {

    std::pair<Real,Real> interval(-1.0, 1.0);
    {
        const size_t n = 1;
        std::pair<Real,Real> alphabeta(1.0, 1.0);
        Rule rule(n, alphabeta, interval);
        std::array<Real,n> expectedPoints = { 0 };

        std::array<Real,n> expectedWeights = { 1.3333334 };

        auto computedPoints  = rule.getPoints();
        auto computedWeights = rule.getWeights();

        for (std::size_t i = 0; i < expectedPoints.size(); ++i) {
            EXPECT_FLOAT_EQ(expectedPoints [i], computedPoints [i]);
            EXPECT_FLOAT_EQ(expectedWeights[i], computedWeights[i]);
        }
    }

    {
        const size_t n = 2;
        std::pair<Real,Real> alphabeta(1.0, 1.0);
        Rule rule(n, alphabeta, interval);
        std::array<Real,n> expectedPoints = {
                -0.447213595499958,
                 0.447213595499958};

        std::array<Real,n> expectedWeights = {
                0.666666666666667,
                0.666666666666667};

        auto computedPoints  = rule.getPoints();
        auto computedWeights = rule.getWeights();

        for (std::size_t i = 0; i < expectedPoints.size(); ++i) {
            EXPECT_FLOAT_EQ(expectedPoints [i], computedPoints [i]);
            EXPECT_FLOAT_EQ(expectedWeights[i], computedWeights[i]);
        }
    }

    {
        const size_t n = 3;
        std::pair<Real,Real> alphabeta(1.0, 1.0);
        Rule rule(n, alphabeta, interval);
        std::array<Real,n> expectedPoints = {
                -0.654653670707977    ,
                 6.982261990806649e-17,
                 0.654653670707977    };

        std::array<Real,n> expectedWeights = {
                0.3111111111111109,
                0.7111111111111111,
                0.3111111111111108};

        auto computedPoints  = rule.getPoints();
        auto computedWeights = rule.getWeights();

        for (std::size_t i = 0; i < expectedPoints.size(); ++i) {
            EXPECT_FLOAT_EQ(expectedPoints [i], computedPoints [i]);
            EXPECT_FLOAT_EQ(expectedWeights[i], computedWeights[i]);
        }
    }

    {
        const size_t n = 8;
        std::pair<Real,Real> alphabeta(1.0, 1.0);

        std::array<Real,n> expectedPoints = {
                -0.9195339081664586,
                -0.7387738651055049,
                -0.4779249498104446,
                -0.1652789576663870,
                 0.1652789576663869,
                 0.4779249498104445,
                 0.7387738651055049,
                 0.9195339081664586};

        std::array<Real,n> expectedWeights = {
                0.02059009564912185,
                0.10214770236035840,
                0.22533655496985810,
                0.31859231368732870,
                0.31859231368732860,
                0.22533655496985820,
                0.10214770236035860,
                0.02059009564912179};

        Rule rule(n, alphabeta, interval);
        auto computedPoints  = rule.getPoints();
        auto computedWeights = rule.getWeights();

        for (std::size_t i = 0; i < expectedPoints.size(); ++i) {
            EXPECT_FLOAT_EQ(expectedPoints [i], computedPoints [i]);
            EXPECT_FLOAT_EQ(expectedWeights[i], computedWeights[i]);
        }
    }

    {
        const size_t n = 8;
        std::pair<Real,Real> alphabeta(0.0, 0.0);

        std::array<Real,n> expectedPoints = {
                -0.9602898564975365,
                -0.7966664774136270,
                -0.5255324099163290,
                -0.1834346424956498,
                 0.1834346424956496,
                 0.5255324099163292,
                 0.7966664774136268,
                 0.9602898564975364};

        std::array<Real,n> expectedWeights = {
                 0.1012285362903760,
                 0.2223810344533743,
                 0.3137066458778873,
                 0.3626837833783622,
                 0.3626837833783619,
                 0.3137066458778869,
                 0.2223810344533742,
                 0.1012285362903759};

        Rule rule(n, alphabeta, interval);
        auto computedPoints  = rule.getPoints();
        auto computedWeights = rule.getWeights();

        for (std::size_t i = 0; i < expectedPoints.size(); ++i) {
            EXPECT_FLOAT_EQ(expectedPoints [i], computedPoints [i]);
            EXPECT_FLOAT_EQ(expectedWeights[i], computedWeights[i]);
        }
    }
}
