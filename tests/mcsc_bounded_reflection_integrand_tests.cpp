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

#include "../src/methods/MCSCBoundedReflectionIntegrand.h"

#include "../src/domains/MCSCCircleDomain.h"
#include "../src/domains/MCSCPolygonalDomain.h"

#include <gtest/gtest.h>
#include <cmath>

namespace
{

// Two-circle bounded configuration: outer circle 0 (unit circle, 3 prevertices -- a degenerate
// "triangle" target so the angle-sum check in MCSCPolygonalDomain passes) and one interior hole
// circle, well inside circle 0. Not a physically realizable MCSC map (no attempt to solve the
// parameter problem here); this is purely a fixture for exercising evalFPrime's arithmetic.
MCSCCircleDomain makeTwoCircleDomain()
{
    return MCSCCircleDomain(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
        {Complex(0.3, 0.0), 0.2, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
    });
}

MCSCPolygonalDomain makeTwoComponentPolygon()
{
    // Two triangles (any simple closed 3-vertex polygon works for the angle-sum check); vertex
    // counts (3 each) must match the circle domain's prevertex counts. The outer component's
    // literal shape doesn't need to enclose the inner one -- MCSCPolygonalDomain validates each
    // component independently.
    std::vector<Complex> outerTriangle = {Complex(0.0, 0.0), Complex(5.0, 0.0), Complex(2.5, 5.0)};
    std::vector<Complex> innerTriangle = {Complex(2.0, 2.0), Complex(3.0, 2.0), Complex(2.5, 3.0)};
    return MCSCPolygonalDomain({outerTriangle, innerTriangle}, /*isUnboundedDomain=*/false);
}

} // namespace

TEST(MCSCBoundedReflectionIntegrandTest, ConstructionSucceedsOnMatchingDomains)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();

    EXPECT_NO_THROW(MCSCBoundedReflectionIntegrand(polygon, circle, 1));
}

TEST(MCSCBoundedReflectionIntegrandTest, ConstructionThrowsOnConnectivityMismatch)
{
    auto circle = makeTwoCircleDomain();
    MCSCPolygonalDomain onePolygon(
        {{Complex(0.0, 0.0), Complex(5.0, 0.0), Complex(2.5, 5.0)}}, /*isUnboundedDomain=*/false);

    EXPECT_THROW(MCSCBoundedReflectionIntegrand(onePolygon, circle, 1), std::invalid_argument);
}

TEST(MCSCBoundedReflectionIntegrandTest, ConstructionThrowsOnVertexCountMismatchOnFirstComponent)
{
    auto circle = makeTwoCircleDomain();  // 3 prevertices per circle
    std::vector<Complex> squareOuter = {
        Complex(0.0, 0.0), Complex(5.0, 0.0), Complex(5.0, 5.0), Complex(0.0, 5.0)};
    std::vector<Complex> innerTriangle = {Complex(2.0, 2.0), Complex(3.0, 2.0), Complex(2.5, 3.0)};
    MCSCPolygonalDomain mismatched({squareOuter, innerTriangle}, /*isUnboundedDomain=*/false);

    EXPECT_THROW(MCSCBoundedReflectionIntegrand(mismatched, circle, 1), std::invalid_argument);
}

TEST(MCSCBoundedReflectionIntegrandTest, ConstructionThrowsOnVertexCountMismatchOnLaterComponent)
{
    // Same as ConstructionThrowsOnVertexCountMismatchOnFirstComponent, but with the mismatch on
    // component 1 instead of component 0 -- validateDomains loops over every component, so a
    // future off-by-one (e.g. stopping the loop one component early) would only be caught by
    // exercising a mismatch that isn't on the first component.
    auto circle = makeTwoCircleDomain();  // 3 prevertices per circle
    std::vector<Complex> outerTriangle = {Complex(0.0, 0.0), Complex(5.0, 0.0), Complex(2.5, 5.0)};
    std::vector<Complex> squareInner = {
        Complex(2.0, 2.0), Complex(3.0, 2.0), Complex(3.0, 3.0), Complex(2.0, 3.0)};
    MCSCPolygonalDomain mismatched({outerTriangle, squareInner}, /*isUnboundedDomain=*/false);

    EXPECT_THROW(MCSCBoundedReflectionIntegrand(mismatched, circle, 1), std::invalid_argument);
}

TEST(MCSCBoundedReflectionIntegrandTest, ConstructionThrowsOnFewerThanTwoCircles)
{
    MCSCCircleDomain oneCircle(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
    });
    MCSCPolygonalDomain onePolygon(
        {{Complex(0.0, 0.0), Complex(5.0, 0.0), Complex(2.5, 5.0)}}, /*isUnboundedDomain=*/false);

    EXPECT_THROW(MCSCBoundedReflectionIntegrand(onePolygon, oneCircle, 1), std::invalid_argument);
}

TEST(MCSCBoundedReflectionIntegrandTest, EvalFPrimeIsFiniteAwayFromPrevertices)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCBoundedReflectionIntegrand integrand(polygon, circle, 2);

    Complex z(0.5, 0.1);  // inside circle 0, away from circle 1 and all reflections
    Complex fp = integrand.evalFPrime(z);

    EXPECT_TRUE(std::isfinite(std::real(fp)));
    EXPECT_TRUE(std::isfinite(std::imag(fp)));
    EXPECT_GT(std::abs(fp), 0.0);
}

TEST(MCSCBoundedReflectionIntegrandTest, EvalFPrimeAtZeroTruncationMatchesHandComputedValue)
{
    // At N=0, evalFPrime's level loop runs zero times (this class's port of fpintrefl.m's
    // `for level=1:N` loop), leaving only the level-0 C1 term -- interior circles contribute
    // nothing at all at N=0 (unlike the unbounded case, where every circle contributes its own
    // level-0 term). This is a strong, distinctive pin: an implementation that mistakenly
    // included interior circles' own level-0 entries at N=0 would fail this test but might still
    // pass a weaker "finite and nonzero" check.
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCBoundedReflectionIntegrand integrand(polygon, circle, 0);

    Complex z(0.5, 0.1);

    const Complex c1(0.0, 0.0);
    std::vector<Complex> c1Prevertices = circle.getPrevertices(0);
    const auto& alpha0 = polygon.getAlpha(0);

    Complex expectedLogSum{0.0, 0.0};
    for (std::size_t k = 0; k < 3; ++k)
    {
        const double beta = -(1.0 - alpha0[k]);
        expectedLogSum += beta * std::log(1.0 - (c1Prevertices[k] - c1) / (z - c1));
    }
    Complex expected = std::exp(expectedLogSum);

    Complex actual = integrand.evalFPrime(z);
    EXPECT_NEAR(std::real(actual), std::real(expected), 1e-10);
    EXPECT_NEAR(std::imag(actual), std::imag(expected), 1e-10);
}

TEST(MCSCBoundedReflectionIntegrandTest, EvalFPrimeVanishesApproachingAC1PrevertexAtExpectedRate)
{
    // Each C1 vertex factor in eval_fprime's product is ((z - z_prevertex)/(z - c1))^beta with
    // beta = -(1 - alpha) for the outer component. For a triangle-ish vertex (0 < alpha < 1, so
    // -1 < beta < 0), this is a pole (not a branch-point zero) -- |f'| ~ C * d^beta -> infinity as
    // z -> z_prevertex. Checking the decay/growth *rate* (not just "large near the prevertex")
    // catches a sign/exponent bug -- e.g. beta accidentally computed as +(1-alpha) instead of
    // -(1-alpha) -- that a simple "large near vertex" check would miss (a branch-point zero at
    // the wrong sign would still show *some* change in magnitude nearby).
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCBoundedReflectionIntegrand integrand(polygon, circle, 1);

    const double beta = -(1.0 - polygon.getAlpha(0)[0]);
    ASSERT_LT(beta, 0.0);
    ASSERT_GT(beta, -1.0);

    Complex prevertex(1.0, 0.0);  // circle 0's prevertex 0, at angle 0 on the unit circle
    const double d1 = 1e-4;
    const double d2 = 1e-2;  // 100x farther
    double mag1 = std::abs(integrand.evalFPrime(prevertex - Complex(d1, 0.0)));
    double mag2 = std::abs(integrand.evalFPrime(prevertex - Complex(d2, 0.0)));

    double observedExponent = std::log(mag2 / mag1) / std::log(d2 / d1);
    EXPECT_NEAR(observedExponent, beta, 0.05);
}

TEST(MCSCBoundedReflectionIntegrandTest, EvalFPrimeVanishesApproachingAnInteriorPrevertexAtExpectedRate)
{
    // Same decay-rate check as above, but for an interior circle's own prevertex -- exercises the
    // ordinary-term branch of the (level, j, nu) loop (and its paired "extra" C1 term) separately
    // from the outer component's level-0 term.
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCBoundedReflectionIntegrand integrand(polygon, circle, 1);

    const double beta = 1.0 - polygon.getAlpha(1)[0];
    ASSERT_GT(beta, 0.0);
    ASSERT_LT(beta, 1.0);

    // Interior circle 1's prevertex 0 is at (0.5, 0.0) (center (0.3, 0.0), radius 0.2, angle 0),
    // but evalFPrime's ordinary term for circle 1 uses its *pre-reflected-to-exterior* copy
    // (mcsc::reflectPoint through C1, the unit circle): c + r^2/conj(z-c) with c=0, r=1 sends
    // (0.5, 0.0) to (2.0, 0.0) -- the singularity this factor actually produces sits there, not
    // at the original prevertex point.
    Complex prevertex(2.0, 0.0);
    const double d1 = 1e-5;
    const double d2 = 1e-3;  // 100x farther
    double mag1 = std::abs(integrand.evalFPrime(prevertex - Complex(d1, 0.0)));
    double mag2 = std::abs(integrand.evalFPrime(prevertex - Complex(d2, 0.0)));

    double observedExponent = std::log(mag2 / mag1) / std::log(d2 / d1);
    EXPECT_NEAR(observedExponent, beta, 0.05);
}

TEST(MCSCBoundedReflectionIntegrandTest, RebuildMatchesFreshConstructionOnUpdatedCircleDomain)
{
    auto circle = makeTwoCircleDomain();
    auto polygon = makeTwoComponentPolygon();
    MCSCBoundedReflectionIntegrand integrand(polygon, circle, 1);

    Complex z(0.5, 0.1);
    Complex before = integrand.evalFPrime(z);

    // Move the interior circle -- rebuild must actually change the result, not silently keep
    // stale reflection data (the whole reason MCSCBoundedReflectionIntegrand requires an explicit
    // rebuild() instead of MATLAB's listener-based auto-rebuild).
    circle.setCenter(1, Complex(-0.2, 0.1));
    integrand.rebuild(polygon, circle, 1);
    Complex after = integrand.evalFPrime(z);
    EXPECT_GT(std::abs(before - after), 1e-9);

    // Not just "changed" -- verify rebuild() produces the *correct* new value by comparing
    // against a fresh MCSCBoundedReflectionIntegrand constructed directly on the moved domain.
    MCSCBoundedReflectionIntegrand freshOnMovedDomain(polygon, circle, 1);
    Complex fresh = freshOnMovedDomain.evalFPrime(z);
    EXPECT_NEAR(std::real(after), std::real(fresh), 1e-12);
    EXPECT_NEAR(std::imag(after), std::imag(fresh), 1e-12);
}
