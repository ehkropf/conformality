/*
 * Copyright © 2026, Everett Kropf (ehkropf@gmail.com)
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

#include "../src/methods/MCSCReflectionIntegrand.h"

#include "../src/domains/MCSCCircleDomain.h"
#include "../src/domains/MCSCPolygonalDomain.h"

#include <gtest/gtest.h>
#include <cmath>

namespace
{

// Two-circle unbounded configuration: outer circle 0 (unit circle, 3 prevertices -- a
// degenerate "triangle" target so the angle-sum check in MCSCPolygonalDomain passes) and one
// hole circle. Not a physically realizable MCSC map (no attempt to solve the parameter problem
// here); this is purely a fixture for exercising evalFPrime/evalFkj's arithmetic.
MCSCCircleDomain makeTwoCircleDomain()
{
    return MCSCCircleDomain(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
        {Complex(3.0, 0.0), 0.5, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
    });
}

MCSCPolygonalDomain makeTwoComponentPolygon()
{
    // Two triangles (any simple closed 3-vertex polygon works for the angle-sum check); vertex
    // counts (3 each) must match the circle domain's prevertex counts.
    std::vector<Complex> triangleA = {Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.5, 1.0)};
    std::vector<Complex> triangleB = {Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(5.5, 1.0)};
    return MCSCPolygonalDomain({triangleA, triangleB}, /*isUnboundedDomain=*/true);
}

} // namespace

TEST(MCSCReflectionIntegrandTest, ConstructionSucceedsOnMatchingDomains)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();

    EXPECT_NO_THROW(MCSCReflectionIntegrand(polygon, circle, 1));
}

TEST(MCSCReflectionIntegrandTest, ConstructionThrowsOnConnectivityMismatch)
{
    auto circle = makeTwoCircleDomain();
    MCSCPolygonalDomain onePolygon(
        {{Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.5, 1.0)}}, /*isUnboundedDomain=*/true);

    EXPECT_THROW(MCSCReflectionIntegrand(onePolygon, circle, 1), std::invalid_argument);
}

TEST(MCSCReflectionIntegrandTest, ConstructionThrowsOnVertexCountMismatch)
{
    auto circle = makeTwoCircleDomain();  // 3 prevertices per circle
    std::vector<Complex> squareA = {
        Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 1.0), Complex(0.0, 1.0)};
    std::vector<Complex> triangleB = {Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(5.5, 1.0)};
    MCSCPolygonalDomain mismatched({squareA, triangleB}, /*isUnboundedDomain=*/true);

    EXPECT_THROW(MCSCReflectionIntegrand(mismatched, circle, 1), std::invalid_argument);
}

TEST(MCSCReflectionIntegrandTest, EvalFPrimeIsFiniteAwayFromPrevertices)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 2);

    Complex z(10.0, 5.0);  // far from both circles and their reflections
    Complex fp = integrand.evalFPrime(z);

    EXPECT_TRUE(std::isfinite(std::real(fp)));
    EXPECT_TRUE(std::isfinite(std::imag(fp)));
    EXPECT_GT(std::abs(fp), 0.0);
}

TEST(MCSCReflectionIntegrandTest, EvalFPrimeVanishesApproachingAPrevertex)
{
    // Each vertex factor in eval_fprime's product is ((z - z_prevertex)/(z - c))^beta with
    // beta = 1 - alpha. For a triangle-ish vertex (0 < alpha < 1, so 0 < beta < 1), this factor
    // -> 0 as z -> z_prevertex (not a pole -- f'(z) itself is proportional to
    // (z - z_prevertex)^beta there, a branch-point zero, matching the standard
    // Schwarz-Christoffel prevertex behavior for f'). Confirms the port carries the right sign
    // convention, not just "doesn't crash."
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 1);

    // Prevertex 0 on circle 0 is at angle 0 on the unit circle, i.e. z = 1 + 0i.
    Complex prevertex(1.0, 0.0);
    Complex near = prevertex + Complex(1e-6, 0.0);
    Complex farther = prevertex + Complex(1e-2, 0.0);

    double magNear = std::abs(integrand.evalFPrime(near));
    double magFarther = std::abs(integrand.evalFPrime(farther));

    EXPECT_LT(magNear, magFarther);
    EXPECT_NEAR(magNear, 0.0, 1e-2);
}

TEST(MCSCReflectionIntegrandTest, EvalFkjThrowsOnOutOfRangeIndices)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 1);

    Complex z(10.0, 5.0);
    EXPECT_THROW(integrand.evalFkj(z, 0, -1), std::invalid_argument);
    EXPECT_THROW(integrand.evalFkj(z, 0, 2), std::invalid_argument);
    EXPECT_THROW(integrand.evalFkj(z, -1, 0), std::invalid_argument);
    EXPECT_THROW(integrand.evalFkj(z, 3, 0), std::invalid_argument);
}

TEST(MCSCReflectionIntegrandTest, EvalFkjIsFiniteAwayFromPrevertices)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 1);

    Complex z(10.0, 5.0);
    Complex fkj = integrand.evalFkj(z, 0, 0);

    EXPECT_TRUE(std::isfinite(std::real(fkj)));
    EXPECT_TRUE(std::isfinite(std::imag(fkj)));
}

TEST(MCSCReflectionIntegrandTest, RebuildReflectsUpdatedCircleDomain)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 1);

    Complex z(10.0, 5.0);
    Complex before = integrand.evalFPrime(z);

    // Move the hole circle -- rebuild must actually change the result, not silently keep stale
    // reflection data (the whole reason MCSCReflectionIntegrand requires an explicit rebuild()
    // instead of MATLAB's listener-based auto-rebuild).
    circle.setCenter(1, Complex(2.0, 1.0));
    integrand.rebuild(polygon, circle, 1);
    Complex after = integrand.evalFPrime(z);

    EXPECT_GT(std::abs(before - after), 1e-9);
}
