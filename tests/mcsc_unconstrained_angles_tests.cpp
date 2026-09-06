/*
 * Copyright © 2025, Everett Kropf (ehkropf@gmail.com)
 *
 * This file is part of Conformality.
 * Conformality is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Conformality is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License along
 * with Conformality. If not, see <https://www.gnu.org/licenses/>.
 */

#include "../src/domains/MCSCUnconstrainedAngles.h"

#include <gtest/gtest.h>
#include <cmath>

using mcsc::anglesFromUnconstrained;
using mcsc::anglesToUnconstrained;

namespace
{
constexpr double TWO_PI = 2.0 * M_PI;
}

TEST(MCSCUnconstrainedAnglesTest, RoundTripPreservesAngles)
{
    std::vector<double> theta = {0.0, 1.0, 2.5, 4.0};

    auto psi = anglesToUnconstrained(theta);
    ASSERT_EQ(psi.size(), theta.size() - 1);

    auto recovered = anglesFromUnconstrained(theta.front(), psi);
    ASSERT_EQ(recovered.size(), theta.size());
    for (std::size_t k = 0; k < theta.size(); ++k)
    {
        EXPECT_NEAR(theta[k], recovered[k], 1e-12);
    }
}

TEST(MCSCUnconstrainedAnglesTest, RoundTripWithNonzeroFirstAngle)
{
    // Mirrors circles j >= 2 in circdomain.m, where theta_1,j is itself a free parameter.
    std::vector<double> theta = {0.7, 1.9, 3.0, 5.5};

    auto psi = anglesToUnconstrained(theta);
    auto recovered = anglesFromUnconstrained(theta.front(), psi);

    ASSERT_EQ(recovered.size(), theta.size());
    for (std::size_t k = 0; k < theta.size(); ++k)
    {
        EXPECT_NEAR(theta[k], recovered[k], 1e-12);
    }
}

TEST(MCSCUnconstrainedAnglesTest, RoundTripMinimalTwoPrevertices)
{
    std::vector<double> theta = {0.0, M_PI};

    auto psi = anglesToUnconstrained(theta);
    ASSERT_EQ(psi.size(), 1u);

    auto recovered = anglesFromUnconstrained(theta.front(), psi);
    ASSERT_EQ(recovered.size(), 2u);
    EXPECT_NEAR(theta[0], recovered[0], 1e-12);
    EXPECT_NEAR(theta[1], recovered[1], 1e-12);
}

TEST(MCSCUnconstrainedAnglesTest, RecoveredAnglesAreAlwaysOrdered)
{
    // Any unconstrained psi vector -- including large-magnitude / mixed-sign entries that would
    // correspond to very unequal gaps -- must still produce strictly increasing angles. This is
    // the entire point of the transform: the ordering constraint is automatically satisfied.
    std::vector<double> psi = {5.0, -5.0, 2.0, -0.3};
    double theta1 = 0.2;

    auto theta = anglesFromUnconstrained(theta1, psi);

    ASSERT_EQ(theta.size(), psi.size() + 1);
    for (std::size_t k = 0; k + 1 < theta.size(); ++k)
    {
        EXPECT_LT(theta[k], theta[k + 1]);
    }
    EXPECT_LT(theta.back(), theta1 + TWO_PI);
    EXPECT_GE(theta.front(), theta1);
}

TEST(MCSCUnconstrainedAnglesTest, GapsSumToTwoPi)
{
    double theta1 = 0.0;
    std::vector<double> psi = {1.3, -0.4, 0.9};

    auto theta = anglesFromUnconstrained(theta1, psi);

    double sum = 0.0;
    for (std::size_t k = 0; k + 1 < theta.size(); ++k)
    {
        sum += theta[k + 1] - theta[k];
    }
    sum += (theta1 + TWO_PI) - theta.back();

    EXPECT_NEAR(sum, TWO_PI, 1e-12);
}

TEST(MCSCUnconstrainedAnglesTest, ThrowsOnTooFewAngles)
{
    std::vector<double> theta = {0.0};
    EXPECT_THROW(anglesToUnconstrained(theta), std::invalid_argument);
}

TEST(MCSCUnconstrainedAnglesTest, ThrowsOnNonIncreasingAngles)
{
    std::vector<double> theta = {0.0, 2.0, 1.0, 4.0};
    EXPECT_THROW(anglesToUnconstrained(theta), std::invalid_argument);
}

TEST(MCSCUnconstrainedAnglesTest, ThrowsWhenAnglesSpanAtLeastTwoPi)
{
    std::vector<double> theta = {0.0, 2.0, 4.0, TWO_PI + 0.1};
    EXPECT_THROW(anglesToUnconstrained(theta), std::invalid_argument);
}
