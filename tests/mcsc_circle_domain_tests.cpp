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

#include "../src/domains/MCSCCircleDomain.h"

#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace
{

std::vector<MCSCCircleDomain::CircleData> makeThreeConnectedExample()
{
    return {
        {Complex(0.0, 0.0), 1.0, {0.0, 1.0, 2.5, 4.0}},
        {Complex(0.2, 0.1), 0.3, {0.5, 2.0, 3.5}},
        {Complex(-0.3, 0.2), 0.2, {1.0, 3.0}},
    };
}

} // namespace

TEST(MCSCCircleDomainTest, ConstructionAndAccessors)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());

    EXPECT_EQ(domain.getConnectivity(), 3);
    EXPECT_EQ(domain.circleCount(), 3);
    EXPECT_FALSE(domain.isUnbounded());

    EXPECT_EQ(domain.getCenter(0), Complex(0.0, 0.0));
    EXPECT_DOUBLE_EQ(domain.getRadius(0), 1.0);
    EXPECT_EQ(domain.prevertexCount(0), 4u);

    EXPECT_EQ(domain.getCenter(1), Complex(0.2, 0.1));
    EXPECT_DOUBLE_EQ(domain.getRadius(1), 0.3);
    EXPECT_EQ(domain.prevertexCount(1), 3u);

    auto centers = domain.getCenters();
    auto radii = domain.getRadii();
    ASSERT_EQ(centers.size(), 3u);
    ASSERT_EQ(radii.size(), 3u);
    EXPECT_EQ(centers[2], Complex(-0.3, 0.2));
    EXPECT_DOUBLE_EQ(radii[2], 0.2);
}

TEST(MCSCCircleDomainTest, ConstructionThrowsOnEmptyCircleList)
{
    EXPECT_THROW(MCSCCircleDomain(std::vector<MCSCCircleDomain::CircleData>{}), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, ConstructionThrowsOnNonPositiveRadius)
{
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 0.0, {0.0, 1.0, 2.0}},
    };
    EXPECT_THROW(MCSCCircleDomain{circles}, std::invalid_argument);
}

TEST(MCSCCircleDomainTest, ConstructionThrowsOnTooFewPrevertices)
{
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0}},
    };
    EXPECT_THROW(MCSCCircleDomain{circles}, std::invalid_argument);
}

TEST(MCSCCircleDomainTest, ConstructionThrowsOnUnorderedPrevertices)
{
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, 2.0, 1.0}},
    };
    EXPECT_THROW(MCSCCircleDomain{circles}, std::invalid_argument);
}

TEST(MCSCCircleDomainTest, GetPrevertices)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());

    auto prevertices = domain.getPrevertices(1);
    ASSERT_EQ(prevertices.size(), 3u);

    const Complex center(0.2, 0.1);
    const double radius = 0.3;
    const std::vector<double> angles = {0.5, 2.0, 3.5};
    for (std::size_t k = 0; k < angles.size(); ++k)
    {
        Complex expected = center + radius * Complex(std::cos(angles[k]), std::sin(angles[k]));
        EXPECT_NEAR(prevertices[k].real(), expected.real(), 1e-14);
        EXPECT_NEAR(prevertices[k].imag(), expected.imag(), 1e-14);
    }
}

TEST(MCSCCircleDomainTest, SettersUpdateState)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());

    domain.setCenter(1, Complex(0.5, -0.1));
    EXPECT_EQ(domain.getCenter(1), Complex(0.5, -0.1));

    domain.setRadius(1, 0.25);
    EXPECT_DOUBLE_EQ(domain.getRadius(1), 0.25);

    domain.setPrevertexAngles(1, {0.1, 1.0, 2.0});
    auto angles = domain.getPrevertexAngles(1);
    ASSERT_EQ(angles.size(), 3u);
    EXPECT_DOUBLE_EQ(angles[0], 0.1);
    EXPECT_DOUBLE_EQ(angles[1], 1.0);
    EXPECT_DOUBLE_EQ(angles[2], 2.0);
}

TEST(MCSCCircleDomainTest, SetRadiusThrowsOnNonPositive)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_THROW(domain.setRadius(1, 0.0), std::invalid_argument);
    EXPECT_THROW(domain.setRadius(1, -1.0), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, SetCenterAndSetRadiusThrowOnCircleZero)
{
    // Circle 0's center/radius are fixed by convention (c_0 = 0, r_0 = 1) -- this convention is
    // relied on by toUnconstrained()/setFromUnconstrained(), which never pack/restore circle 0's
    // center or radius, so mutating them would otherwise be silently undone on the next round trip.
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_THROW(domain.setCenter(0, Complex(0.1, 0.1)), std::invalid_argument);
    EXPECT_THROW(domain.setRadius(0, 0.5), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, SetPrevertexAnglesThrowsOnCountMismatch)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_THROW(domain.setPrevertexAngles(0, {0.0, 1.0}), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, SetPrevertexAnglesThrowsOnUnordered)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_THROW(domain.setPrevertexAngles(2, {2.0, 0.5}), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, IndexAccessorsThrowOutOfRange)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_THROW(domain.getCenter(-1), std::invalid_argument);
    EXPECT_THROW(domain.getCenter(3), std::invalid_argument);
    EXPECT_THROW(domain.getRadius(3), std::invalid_argument);
    EXPECT_THROW(domain.getPrevertexAngles(3), std::invalid_argument);
    EXPECT_THROW(domain.getPrevertices(3), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, ContainsInsideOuterOutsideHoles)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());

    // (0.9, 0.0) is inside the outer unit circle and outside both holes (holes are centered near
    // the origin with small radii, so the origin itself lies inside hole 1 -- avoid it here).
    EXPECT_TRUE(domain.contains(Complex(0.9, 0.0)));

    // Center of a hole is not in the domain.
    EXPECT_FALSE(domain.contains(Complex(0.2, 0.1)));
    EXPECT_FALSE(domain.contains(Complex(-0.3, 0.2)));

    // Outside the outer circle entirely.
    EXPECT_FALSE(domain.contains(Complex(2.0, 0.0)));
}

TEST(MCSCCircleDomainTest, ContainsIsFalseExactlyOnBoundary)
{
    // Open-domain semantics, no boundary tolerance: a point exactly on the outer circle is not
    // contained (see contains()'s docstring).
    MCSCCircleDomain domain(makeThreeConnectedExample());
    EXPECT_FALSE(domain.contains(Complex(1.0, 0.0)));
}

TEST(MCSCCircleDomainTest, UnconstrainedRoundTrip)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    MCSCCircleDomain reference(makeThreeConnectedExample());

    auto Xu = domain.toUnconstrained();
    ASSERT_EQ(Xu.size(), 3 * 3 - 4 + 4 + 3 + 2);

    // Round trip through a domain that starts from the same shape but is otherwise independent.
    reference.setFromUnconstrained(Xu);

    for (int j = 0; j < domain.circleCount(); ++j)
    {
        EXPECT_NEAR(reference.getCenter(j).real(), domain.getCenter(j).real(), 1e-12);
        EXPECT_NEAR(reference.getCenter(j).imag(), domain.getCenter(j).imag(), 1e-12);
        EXPECT_NEAR(reference.getRadius(j), domain.getRadius(j), 1e-12);

        auto expectedAngles = domain.getPrevertexAngles(j);
        auto actualAngles = reference.getPrevertexAngles(j);
        ASSERT_EQ(actualAngles.size(), expectedAngles.size());
        for (std::size_t k = 0; k < expectedAngles.size(); ++k)
        {
            EXPECT_NEAR(actualAngles[k], expectedAngles[k], 1e-12);
        }
    }
}

TEST(MCSCCircleDomainTest, UnconstrainedRoundTripAfterMutation)
{
    // Pack -> unpack should also work as a no-op fixed point when applied to the domain's own
    // current state after mutating it in place (the actual usage pattern during a Newton solve).
    MCSCCircleDomain domain(makeThreeConnectedExample());

    domain.setCenter(1, Complex(0.4, -0.15));
    domain.setRadius(2, 0.15);
    domain.setPrevertexAngles(0, {0.2, 1.5, 3.0, 5.0});

    auto Xu = domain.toUnconstrained();
    domain.setFromUnconstrained(Xu);
    auto Xu2 = domain.toUnconstrained();

    ASSERT_EQ(Xu.size(), Xu2.size());
    for (int i = 0; i < Xu.size(); ++i)
    {
        EXPECT_NEAR(Xu[i], Xu2[i], 1e-10);
    }
}

TEST(MCSCCircleDomainTest, SetFromUnconstrainedThrowsOnWrongSize)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    Eigen::VectorXd tooShort = Eigen::VectorXd::Zero(3);
    EXPECT_THROW(domain.setFromUnconstrained(tooShort), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, TransformBoundaryTranslatesCenters)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    Complex shift(1.0, -2.0);

    domain.transformBoundary([shift](const Complex& z) { return z + shift; });

    EXPECT_NEAR(domain.getCenter(0).real(), shift.real(), 1e-12);
    EXPECT_NEAR(domain.getCenter(0).imag(), shift.imag(), 1e-12);
    EXPECT_NEAR(domain.getRadius(0), 1.0, 1e-12);

    EXPECT_NEAR(domain.getCenter(1).real(), 0.2 + shift.real(), 1e-12);
    EXPECT_NEAR(domain.getCenter(1).imag(), 0.1 + shift.imag(), 1e-12);
    EXPECT_NEAR(domain.getRadius(1), 0.3, 1e-12);
}

TEST(MCSCCircleDomainTest, TransformBoundaryScalesRadii)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());

    domain.transformBoundary([](const Complex& z) { return 2.0 * z; });

    EXPECT_NEAR(domain.getRadius(0), 2.0, 1e-12);
    EXPECT_NEAR(domain.getRadius(1), 0.6, 1e-12);
    EXPECT_NEAR(domain.getCenter(1).real(), 0.4, 1e-12);
    EXPECT_NEAR(domain.getCenter(1).imag(), 0.2, 1e-12);
}

TEST(MCSCCircleDomainTest, ConstructionThrowsWhenPrevertexSpanIsExactlyTwoPi)
{
    // A span of exactly 2*pi is rejected; only strictly less than 2*pi is valid (see
    // validateCircleData()'s t.back() >= t.front() + 2*pi check).
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, M_PI, 2.0 * M_PI}},
    };
    EXPECT_THROW(MCSCCircleDomain{circles}, std::invalid_argument);
}

TEST(MCSCCircleDomainTest, SimplyConnectedSingleCircle)
{
    // m = 1 is the smallest legal circle count per the class's own docstring ("size m >= 1"); the
    // unconstrained vector then contains only circle 0's angle block (no radii/center blocks).
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, 1.0, 2.5, 4.0}},
    };
    MCSCCircleDomain domain(circles);
    MCSCCircleDomain reference(circles);

    EXPECT_EQ(domain.circleCount(), 1);

    auto Xu = domain.toUnconstrained();
    ASSERT_EQ(Xu.size(), 4 - 1); // K0 - 1 unconstrained angle variables, no radii/center blocks.

    reference.setFromUnconstrained(Xu);
    auto expectedAngles = domain.getPrevertexAngles(0);
    auto actualAngles = reference.getPrevertexAngles(0);
    ASSERT_EQ(actualAngles.size(), expectedAngles.size());
    for (std::size_t k = 0; k < expectedAngles.size(); ++k)
    {
        EXPECT_NEAR(actualAngles[k], expectedAngles[k], 1e-12);
    }
}

TEST(MCSCCircleDomainTest, FourConnectedRoundTrip)
{
    // Connectivity 4 exercises the center-block index arithmetic (base = (m-1) + 2*(j-1)) beyond
    // the smallest cases, per CLAUDE.md's guidance to test loop-based initialization at higher
    // connectivity, not just 2 or 3.
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, 1.0, 2.5, 4.0}},
        {Complex(0.3, 0.1), 0.2, {0.5, 2.0, 3.5}},
        {Complex(-0.3, 0.2), 0.15, {1.0, 3.0}},
        {Complex(0.1, -0.4), 0.1, {0.2, 1.8, 3.1, 4.5}},
    };
    MCSCCircleDomain domain(circles);
    MCSCCircleDomain reference(circles);

    auto Xu = domain.toUnconstrained();
    ASSERT_EQ(Xu.size(), 3 * 4 - 4 + 4 + 3 + 2 + 4);

    reference.setFromUnconstrained(Xu);
    for (int j = 0; j < domain.circleCount(); ++j)
    {
        EXPECT_NEAR(reference.getCenter(j).real(), domain.getCenter(j).real(), 1e-12);
        EXPECT_NEAR(reference.getCenter(j).imag(), domain.getCenter(j).imag(), 1e-12);
        EXPECT_NEAR(reference.getRadius(j), domain.getRadius(j), 1e-12);

        auto expectedAngles = domain.getPrevertexAngles(j);
        auto actualAngles = reference.getPrevertexAngles(j);
        ASSERT_EQ(actualAngles.size(), expectedAngles.size());
        for (std::size_t k = 0; k < expectedAngles.size(); ++k)
        {
            EXPECT_NEAR(actualAngles[k], expectedAngles[k], 1e-12);
        }
    }
}

TEST(MCSCCircleDomainTest, UnconstrainedLayoutMatchesDocumentedPackingOrder)
{
    // Pins the raw Xu layout against the packing order documented on toUnconstrained(), so a
    // block reordering that happens to still round-trip internally (pack and unpack agreeing with
    // each other, but not with the documented/MATLAB-matching contract) would be caught here.
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, M_PI}},           // circle 0: K0 = 2, so 1 psi value.
        {Complex(0.5, -0.25), 2.0, {0.25 * M_PI, M_PI}}, // circle 1: K1 = 2, so 1 psi value.
    };
    MCSCCircleDomain domain(circles);

    auto Xu = domain.toUnconstrained();
    ASSERT_EQ(Xu.size(), 6); // 3*2-4 + K0 + K1 = 2 + 2 + 2 = 6.

    // Block 1 (1 entry, m-1=1): log(r_1).
    EXPECT_NEAR(Xu[0], std::log(2.0), 1e-12);
    // Block 2 (2 entries, 2*(m-1)=2): Re(c_1), Im(c_1).
    EXPECT_NEAR(Xu[1], 0.5, 1e-12);
    EXPECT_NEAR(Xu[2], -0.25, 1e-12);
    // Block 3 (K0-1=1 entry): psi for circle 0 (K0 = 2 -> phi_1 = pi, phi_2 = pi,
    // psi_1 = log(phi_2/phi_1) = 0).
    EXPECT_NEAR(Xu[3], 0.0, 1e-12);
    // Block 4 (K1=2 entries): t(1,1), then psi for circle 1 (K1 = 2 -> phi_1 = 0.75*pi,
    // phi_2 = 1.25*pi, psi_1 = log(phi_2/phi_1)).
    EXPECT_NEAR(Xu[4], 0.25 * M_PI, 1e-12);
    EXPECT_NEAR(Xu[5], std::log(1.25 / 0.75), 1e-12);
}

TEST(MCSCCircleDomainTest, SetFromUnconstrainedThrowsOnNonFiniteInput)
{
    MCSCCircleDomain domain(makeThreeConnectedExample());
    auto Xu = domain.toUnconstrained();

    Eigen::VectorXd badXu = Xu;
    badXu[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(domain.setFromUnconstrained(badXu), std::invalid_argument);

    Eigen::VectorXd infXu = Xu;
    infXu[1] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(domain.setFromUnconstrained(infXu), std::invalid_argument);
}

TEST(MCSCCircleDomainTest, TwoConnectedMinimalCase)
{
    // Connectivity 2: no interleaved-center block (m-1 = 1 circle to pack).
    std::vector<MCSCCircleDomain::CircleData> circles = {
        {Complex(0.0, 0.0), 1.0, {0.0, M_PI}},
        {Complex(0.1, -0.1), 0.4, {0.3, 2.0, 4.0}},
    };
    MCSCCircleDomain domain(circles);
    MCSCCircleDomain reference(circles);

    auto Xu = domain.toUnconstrained();
    ASSERT_EQ(Xu.size(), 3 * 2 - 4 + 2 + 3);

    reference.setFromUnconstrained(Xu);
    EXPECT_NEAR(reference.getRadius(1), domain.getRadius(1), 1e-12);
    EXPECT_NEAR(reference.getCenter(1).real(), domain.getCenter(1).real(), 1e-12);
    EXPECT_NEAR(reference.getCenter(1).imag(), domain.getCenter(1).imag(), 1e-12);
}
