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

TEST(MCSCReflectionIntegrandTest, ConstructionThrowsOnVertexCountMismatchOnFirstComponent)
{
    auto circle = makeTwoCircleDomain();  // 3 prevertices per circle
    std::vector<Complex> squareA = {
        Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 1.0), Complex(0.0, 1.0)};
    std::vector<Complex> triangleB = {Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(5.5, 1.0)};
    MCSCPolygonalDomain mismatched({squareA, triangleB}, /*isUnboundedDomain=*/true);

    EXPECT_THROW(MCSCReflectionIntegrand(mismatched, circle, 1), std::invalid_argument);
}

TEST(MCSCReflectionIntegrandTest, ConstructionThrowsOnVertexCountMismatchOnLaterComponent)
{
    // Same as ConstructionThrowsOnVertexCountMismatchOnFirstComponent, but with the mismatch on
    // component 1 instead of component 0 -- validateDomains loops over every component, so a
    // future off-by-one (e.g. stopping the loop one component early) would only be caught by
    // exercising a mismatch that isn't on the first component.
    auto circle = makeTwoCircleDomain();  // 3 prevertices per circle
    std::vector<Complex> triangleA = {Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.5, 1.0)};
    std::vector<Complex> squareB = {
        Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(6.0, 1.0), Complex(5.0, 1.0)};
    MCSCPolygonalDomain mismatched({triangleA, squareB}, /*isUnboundedDomain=*/true);

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

TEST(MCSCReflectionIntegrandTest, EvalFPrimeAtZeroTruncationMatchesHandComputedValue)
{
    // At N=0, reflectCircleSequence produces exactly one (unreflected) entry per circle, with
    // outerImage == center for the unbounded case -- so the ((z-c)/(z-outerImage))^2 factor is
    // identically 1 and evalFPrime reduces to a plain product over the two circles' own
    // prevertices: exp(sum_j sum_k beta[j][k] * log(1 - (z_{k,j} - c_j)/(z - c_j))). This is
    // small enough to compute independently here, pinning the base case that N=1/2 tests (which
    // include reflections) cannot isolate -- an off-by-one that accidentally dropped or
    // duplicated the level-0 term would still "look right" under those.
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 0);

    Complex z(10.0, 5.0);

    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0)};
    std::vector<std::vector<Complex>> prevertices = {
        circle.getPrevertices(0),
        circle.getPrevertices(1),
    };
    std::vector<std::vector<double>> alpha = {polygon.getAlpha(0), polygon.getAlpha(1)};

    Complex expectedLogSum{0.0, 0.0};
    for (int j = 0; j < 2; ++j)
    {
        for (std::size_t k = 0; k < 3; ++k)
        {
            const double beta = 1.0 - alpha[j][k];
            expectedLogSum += beta * std::log(1.0 - (prevertices[j][k] - centers[j]) / (z - centers[j]));
        }
    }
    Complex expected = std::exp(expectedLogSum);

    Complex actual = integrand.evalFPrime(z);
    EXPECT_NEAR(std::real(actual), std::real(expected), 1e-10);
    EXPECT_NEAR(std::imag(actual), std::imag(expected), 1e-10);
}

TEST(MCSCReflectionIntegrandTest, EvalFPrimeVanishesApproachingAPrevertexAtExpectedRate)
{
    // Each vertex factor in eval_fprime's product is ((z - z_prevertex)/(z - c))^beta with
    // beta = 1 - alpha. For a triangle-ish vertex (0 < alpha < 1, so 0 < beta < 1), this factor
    // -> 0 as z -> z_prevertex (not a pole -- f'(z) itself is proportional to
    // (z - z_prevertex)^beta there, a branch-point zero, matching the standard
    // Schwarz-Christoffel prevertex behavior for f'). Checking the decay *rate* (not just that
    // some decay happens) catches a sign/exponent bug -- e.g. beta accidentally computed as
    // alpha instead of 1 - alpha -- that a simple "smaller than farther away" check would miss,
    // since many turning angles would still show some decay either way.
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCReflectionIntegrand integrand(polygon, circle, 1);

    // Prevertex 0 on circle 0 is at angle 0 on the unit circle, i.e. z = 1 + 0i. beta for this
    // vertex: alpha is the interior angle of an equilateral-ish triangleA vertex at (0,0);
    // compute it the same way MCSCPolygonalDomain does, via getAlpha(), to avoid hand-deriving
    // the exact turning angle here.
    const double beta = 1.0 - polygon.getAlpha(0)[0];
    ASSERT_GT(beta, 0.0);
    ASSERT_LT(beta, 1.0);

    Complex prevertex(1.0, 0.0);
    const double d1 = 1e-4;
    const double d2 = 1e-2;  // 100x farther
    double mag1 = std::abs(integrand.evalFPrime(prevertex + Complex(d1, 0.0)));
    double mag2 = std::abs(integrand.evalFPrime(prevertex + Complex(d2, 0.0)));

    // |f'| ~ C * d^beta near the prevertex (the other, non-vanishing factors are locally
    // ~constant over this small a neighborhood), so log(mag2/mag1) / log(d2/d1) should
    // approximate beta.
    double observedExponent = std::log(mag2 / mag1) / std::log(d2 / d1);
    EXPECT_NEAR(observedExponent, beta, 0.05);
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

TEST(MCSCReflectionIntegrandTest, RebuildMatchesFreshConstructionOnUpdatedCircleDomain)
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

    // Not just "changed" -- verify rebuild() produces the *correct* new value by comparing
    // against a fresh MCSCReflectionIntegrand constructed directly on the moved domain. A
    // rebuild() that scrambled the reflection data some other, wrong way would still pass the
    // "before != after" check above but fail this one.
    MCSCReflectionIntegrand freshOnMovedDomain(polygon, circle, 1);
    Complex fresh = freshOnMovedDomain.evalFPrime(z);
    EXPECT_NEAR(std::real(after), std::real(fresh), 1e-12);
    EXPECT_NEAR(std::imag(after), std::imag(fresh), 1e-12);
}
