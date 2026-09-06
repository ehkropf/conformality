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

#include "../src/numerics/GaussJacobiQuadrature.h"

#include <gtest/gtest.h>
#include <cmath>

TEST(GaussJacobiQuadratureTest, LegendreNodesForN1)
{
    auto [nodes, weights] = GaussJacobiQuadrature::gaussianJacobiNodes(1, 0.0, 0.0);
    ASSERT_EQ(nodes.size(), 1);
    EXPECT_NEAR(nodes[0], 0.0, 1e-12);
    EXPECT_NEAR(weights[0], 2.0, 1e-12);
}

TEST(GaussJacobiQuadratureTest, LegendreNodesForN2)
{
    auto [nodes, weights] = GaussJacobiQuadrature::gaussianJacobiNodes(2, 0.0, 0.0);
    ASSERT_EQ(nodes.size(), 2);
    const double expected = 1.0 / std::sqrt(3.0);
    EXPECT_NEAR(nodes[0], -expected, 1e-10);
    EXPECT_NEAR(nodes[1], expected, 1e-10);
    EXPECT_NEAR(weights[0], 1.0, 1e-10);
    EXPECT_NEAR(weights[1], 1.0, 1e-10);
}

TEST(GaussJacobiQuadratureTest, NodesThrowOnInvalidN)
{
    EXPECT_THROW(GaussJacobiQuadrature::gaussianJacobiNodes(0, 0.0, 0.0), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, NodesThrowOnInvalidExponents)
{
    EXPECT_THROW(GaussJacobiQuadrature::gaussianJacobiNodes(4, -1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(GaussJacobiQuadrature::gaussianJacobiNodes(4, 0.0, -1.5), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, IntegratesLowDegreePolynomialExactlyForLegendreCase)
{
    // Gauss-Legendre quadrature with n nodes is exact for polynomials up to degree 2n-1.
    // n=4, alf=bet=0: exact for degree <= 7. Integral of x^6 over [-1,1] is 2/7.
    auto [nodes, weights] = GaussJacobiQuadrature::gaussianJacobiNodes(4, 0.0, 0.0);
    double sum = 0.0;
    for (int i = 0; i < nodes.size(); ++i)
    {
        sum += weights[i] * std::pow(nodes[i], 6);
    }
    EXPECT_NEAR(sum, 2.0 / 7.0, 1e-10);
}

TEST(GaussJacobiQuadratureTest, IntegratesExactlyWithNonzeroBetaWeight)
{
    // integral_{-1}^{1} (1+x)^1 dx = 2 (a constant integrand times the weight function itself,
    // exactly representable by a degree-0 polynomial -- exact for any n >= 1).
    auto [nodes, weights] = GaussJacobiQuadrature::gaussianJacobiNodes(4, 0.0, 1.0);
    double sum = 0.0;
    for (int i = 0; i < nodes.size(); ++i)
    {
        sum += weights[i];
    }
    EXPECT_NEAR(sum, 2.0, 1e-10);
}

TEST(GaussJacobiQuadratureTest, ConstructorFloorsNgjToFour)
{
    // ngj=1 should be floored to 4 -- verify indirectly via node count on the ordinary table by
    // integrating a degree-7 polynomial exactly (only possible with >= 4 nodes).
    GaussJacobiQuadrature quad(std::vector<double>{}, /*ngj=*/1);
    Complex result = quad.integrateLine(
        Complex(-1.0, 0.0), -1, Complex(1.0, 0.0), -1,
        [](const Complex& z) { return std::pow(z, 6); }, {});
    EXPECT_NEAR(std::real(result), 2.0 / 7.0, 1e-8);
}

TEST(GaussJacobiQuadratureTest, IntegrateLineOfConstantFunctionGivesSegmentLength)
{
    GaussJacobiQuadrature quad(std::vector<double>{});
    Complex result = quad.integrateLine(
        Complex(1.0, 1.0), -1, Complex(4.0, 5.0), -1, [](const Complex&) { return Complex(1.0, 0.0); }, {});
    EXPECT_NEAR(std::abs(result), 5.0, 1e-9);  // |4+5i - (1+i)| = |3+4i| = 5
}

TEST(GaussJacobiQuadratureTest, IntegrateArcOfConstantFunctionGivesComplexArcIntegral)
{
    // integral over the unit circle arc from angle 0 to pi/2 of f(z)=1 dz
    // = integral_0^{pi/2} i*exp(i*theta) dtheta = [exp(i*theta)]_0^{pi/2} = exp(i*pi/2) - 1
    // = i - 1 = -1 + i.
    GaussJacobiQuadrature quad(std::vector<double>{});
    Complex result = quad.integrateArc(
        0.0, -1, M_PI / 2.0, -1, Complex(0.0, 0.0), 1.0, [](const Complex&) { return Complex(1.0, 0.0); }, {});
    EXPECT_NEAR(std::real(result), -1.0, 1e-8);
    EXPECT_NEAR(std::imag(result), 1.0, 1e-8);
}

namespace
{
Complex constantOne(const Complex&)
{
    return Complex(1.0, 0.0);
}
}  // namespace

TEST(GaussJacobiQuadratureTest, IntegrateLineThrowsOnInvalidVertexIndex)
{
    GaussJacobiQuadrature quad(std::vector<double>{0.5});
    EXPECT_THROW(
        quad.integrateLine(Complex(0.0, 0.0), 5, Complex(1.0, 0.0), -1, constantOne, {}), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, IntegrateLineThrowsWhenVertexBetaIsAtOrBelowMinusOne)
{
    // beta = -1 has no valid Gauss-Jacobi weight function (matches gjquad.m's calc_qdata
    // skipping beta <= -1) -- using vertex 0 as a singular endpoint should throw.
    GaussJacobiQuadrature quad(std::vector<double>{-1.0});
    EXPECT_THROW(
        quad.integrateLine(Complex(0.0, 0.0), 0, Complex(1.0, 0.0), -1, constantOne, {}), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, IntegrateArcThrowsOnNonPositiveRadius)
{
    GaussJacobiQuadrature quad(std::vector<double>{});
    EXPECT_THROW(
        quad.integrateArc(0.0, -1, 1.0, -1, Complex(0.0, 0.0), 0.0, constantOne, {}), std::invalid_argument);
    EXPECT_THROW(
        quad.integrateArc(0.0, -1, 1.0, -1, Complex(0.0, 0.0), -1.0, constantOne, {}), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, IntegratesAlgebraicSingularityAtLeftEndpointExactly)
{
    // integral_0^1 (x)^(-0.5) dx = 2 -- a genuine algebraic singularity at the left endpoint,
    // the actual point of Gauss-Jacobi quadrature (beta = -0.5 for the (1+t) factor after
    // shifting/scaling the [-1,1] Jacobi weight onto [0,1]; here it's simplest to construct the
    // singularity directly at the left endpoint of a real-axis line segment).
    GaussJacobiQuadrature quad(std::vector<double>{-0.5});
    Complex result = quad.integrateLine(
        Complex(0.0, 0.0), 0, Complex(1.0, 0.0), -1,
        [](const Complex& z) { return std::pow(z, -0.5); }, {});
    EXPECT_NEAR(std::real(result), 2.0, 1e-6);
    EXPECT_NEAR(std::imag(result), 0.0, 1e-6);
}

TEST(GaussJacobiQuadratureTest, IntegratesAlgebraicSingularityAtRightEndpointExactly)
{
    // Same integral as above but with the singularity at the *right* endpoint instead --
    // integral_0^1 (1-x)^(-0.5) dx = 2, exercising the rpn-triggered midpoint split.
    GaussJacobiQuadrature quad(std::vector<double>{-0.5});
    Complex result = quad.integrateLine(
        Complex(0.0, 0.0), -1, Complex(1.0, 0.0), 0,
        [](const Complex& z) { return std::pow(1.0 - z, -0.5); }, {});
    EXPECT_NEAR(std::real(result), 2.0, 1e-6);
    EXPECT_NEAR(std::imag(result), 0.0, 1e-6);
}

TEST(GaussJacobiQuadratureTest, OneHalfRuleHandlesNearbySingularityAccurately)
{
    // Integrate f(z) = 1/(z - p) along a segment with a tracked singularity p close to the
    // *left endpoint* -- the one-half rule's step-length check is evaluated incrementally from
    // the current walking point (matches gjquad.m's lineq/arcq: `abs(a - sng)`, not a lookahead
    // over the whole candidate subinterval), so it only reliably catches a singularity that is
    // near wherever the walk currently stands -- i.e. near an endpoint, which is exactly the
    // case that matters for MCSC (nearby reflected prevertices cluster near path endpoints, not
    // path interiors). A singularity near the *middle* of a long segment is NOT reliably caught
    // by this check on the first step (verified separately, not asserted here) -- a known,
    // faithfully-ported limitation of the one-half rule as specified in gjquad.m, not a bug in
    // this port.
    //
    // Integral of 1/(x-p) dx from 0 to 1: known closed form log((1-p)/(-p)) (verified
    // independently against brute-force numerical integration).
    Complex p(0.05, 0.01);
    GaussJacobiQuadrature quad(std::vector<double>{});
    Complex result = quad.integrateLine(
        Complex(0.0, 0.0), -1, Complex(1.0, 0.0), -1,
        [p](const Complex& z) { return 1.0 / (z - p); }, {p});

    Complex expected = std::log((1.0 - p) / (-p));
    EXPECT_NEAR(std::real(result), std::real(expected), 0.05);
    EXPECT_NEAR(std::imag(result), std::imag(expected), 0.05);

    // Sanity check that this test actually exercises subdivision (i.e. isn't accidentally
    // testing the same single-shot behavior as the "no subdivision needed" case above): a
    // plain, non-adaptive 12-point Gauss-Legendre pass over the same integral is measurably
    // less accurate than the adaptive result.
    GaussJacobiQuadrature noHalfRule(std::vector<double>{}, 12, /*useHalfRule=*/false);
    Complex singleShot = noHalfRule.integrateLine(
        Complex(0.0, 0.0), -1, Complex(1.0, 0.0), -1, [p](const Complex& z) { return 1.0 / (z - p); }, {});
    EXPECT_GT(std::abs(singleShot - expected), std::abs(result - expected));
}

TEST(GaussJacobiQuadratureTest, DisablingHalfRuleStillIntegratesSmoothFunctions)
{
    GaussJacobiQuadrature quad(std::vector<double>{}, 12, /*useHalfRule=*/false);
    Complex result = quad.integrateLine(
        Complex(-1.0, 0.0), -1, Complex(1.0, 0.0), -1, [](const Complex& z) { return std::pow(z, 6); }, {});
    EXPECT_NEAR(std::real(result), 2.0 / 7.0, 1e-8);
}

TEST(GaussJacobiQuadratureTest, DisablingHalfRuleStillHandlesEndpointSingularity)
{
    // useHalfRule only gates the step-length computation, not vertex-table selection -- a
    // genuine endpoint singularity should still be handled correctly with the one-half rule
    // off, in a single Gauss-Jacobi step.
    GaussJacobiQuadrature quad(std::vector<double>{-0.5}, 12, /*useHalfRule=*/false);
    Complex result = quad.integrateLine(
        Complex(0.0, 0.0), 0, Complex(1.0, 0.0), -1, [](const Complex& z) { return std::pow(z, -0.5); }, {});
    EXPECT_NEAR(std::real(result), 2.0, 1e-6);
    EXPECT_NEAR(std::imag(result), 0.0, 1e-6);
}

TEST(GaussJacobiQuadratureTest, IntegratesWithBothEndpointsSingularSimultaneously)
{
    // integral_0^1 x^(-0.5) * (1-x)^(-0.5) dx = B(0.5, 0.5) = Gamma(0.5)^2 / Gamma(1) = pi -- a
    // symmetric Beta-function integral with a genuine algebraic singularity at *both*
    // endpoints simultaneously, exercising integrateWithHalfRule's only recursive path (the
    // midpoint-split-and-subtract triggered by rightVertex >= 0), which the earlier
    // single-endpoint tests above never reach.
    GaussJacobiQuadrature quad(std::vector<double>{-0.5, -0.5});
    Complex result = quad.integrateLine(
        Complex(0.0, 0.0), 0, Complex(1.0, 0.0), 1,
        [](const Complex& z) { return std::pow(z, -0.5) * std::pow(1.0 - z, -0.5); }, {});
    EXPECT_NEAR(std::real(result), M_PI, 1e-6);
    EXPECT_NEAR(std::imag(result), 0.0, 1e-6);
}

TEST(GaussJacobiQuadratureTest, IntegrateArcOfConstantFunctionWithNonUnitRadius)
{
    // Same structural integral as IntegrateArcOfConstantFunctionGivesComplexArcIntegral but at
    // r=2.5 instead of r=1 -- exercises the arcLengthScale unification with a non-trivial
    // scale factor, since at r=1 a transposed r/(1/r) bug would be invisible.
    // integral over the arc (center 0, radius r) from angle 0 to pi/2 of f(z)=1 dz
    // = integral_0^{pi/2} i*r*exp(i*theta) dtheta = r*(exp(i*pi/2) - 1) = r*(-1+i).
    const double r = 2.5;
    GaussJacobiQuadrature quad(std::vector<double>{});
    Complex result = quad.integrateArc(
        0.0, -1, M_PI / 2.0, -1, Complex(0.0, 0.0), r, [](const Complex&) { return Complex(1.0, 0.0); }, {});
    EXPECT_NEAR(std::real(result), r * -1.0, 1e-7);
    EXPECT_NEAR(std::imag(result), r * 1.0, 1e-7);
}

TEST(GaussJacobiQuadratureTest, IntegrateArcThrowsOnInvalidVertexIndex)
{
    // integrateArc delegates vertex validation to the same shared helper as integrateLine;
    // exercise it through this entry point too so a future divergence (e.g. a fast path added
    // only to integrateArc) would be caught.
    GaussJacobiQuadrature quad(std::vector<double>{0.5});
    EXPECT_THROW(
        quad.integrateArc(0.0, 5, M_PI / 2.0, -1, Complex(0.0, 0.0), 1.0, constantOne, {}), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, IntegrateArcThrowsWhenVertexBetaIsAtOrBelowMinusOne)
{
    GaussJacobiQuadrature quad(std::vector<double>{-1.0});
    EXPECT_THROW(
        quad.integrateArc(0.0, 0, M_PI / 2.0, -1, Complex(0.0, 0.0), 1.0, constantOne, {}), std::invalid_argument);
}

TEST(GaussJacobiQuadratureTest, ThrowsConvergenceErrorWhenSubdivisionCapExceeded)
{
    // A dense, evenly-spaced cluster of tracked singularities spanning the segment forces every
    // step to be capped at roughly 2x the (small, roughly constant) spacing -- since the
    // one-half rule only ever allows a step up to 2x the distance to the *nearest* singularity,
    // the walk can never take a step larger than about 2x the cluster spacing, requiring far
    // more than 100 steps to cross a unit-length segment with a spacing this fine.
    std::vector<Complex> singularities;
    for (int i = 1; i <= 300; ++i)
    {
        singularities.emplace_back(i * 0.003, 0.0);
    }
    GaussJacobiQuadrature quad(std::vector<double>{});
    EXPECT_THROW(
        quad.integrateLine(Complex(0.0, 0.0), -1, Complex(1.0, 0.0), -1, constantOne, singularities),
        GaussJacobiQuadrature::ConvergenceError);
}

TEST(GaussJacobiQuadratureTest, IntegratesCorrectlyWhenSingularityListIncludesTheStartingPrevertex)
{
    // Regression test: MCSC's real usage integrates a path starting exactly at a prevertex
    // while passing the *entire* prevertex list (including that same starting prevertex) as
    // the tracked singularity set (see fpextrefl.m's I.sing, populated with every prevertex).
    // A singularity exactly coincident with the current walking position must not be treated as
    // "nearby" -- otherwise the reported distance is 0, collapsing every candidate step length
    // to 0 and spinning the subdivision loop to its cap instead of ever making progress.
    GaussJacobiQuadrature quad(std::vector<double>{-0.5});
    Complex result;
    EXPECT_NO_THROW(
        result = quad.integrateLine(
            Complex(0.0, 0.0), 0, Complex(1.0, 0.0), -1, [](const Complex& z) { return std::pow(z, -0.5); },
            {Complex(0.0, 0.0)}));
    EXPECT_NEAR(std::real(result), 2.0, 1e-6);
    EXPECT_NEAR(std::imag(result), 0.0, 1e-6);
}

TEST(GaussJacobiQuadratureTest, IntegratesArcCorrectlyWhenSingularityListIncludesTheStartingPrevertex)
{
    // Same regression as above, for integrateArc: starting at angle 0 on the unit circle (point
    // 1+0i) with that same point included in the tracked singularity list.
    GaussJacobiQuadrature quad(std::vector<double>{});
    Complex result;
    EXPECT_NO_THROW(
        result = quad.integrateArc(
            0.0, -1, M_PI / 2.0, -1, Complex(0.0, 0.0), 1.0, [](const Complex&) { return Complex(1.0, 0.0); },
            {Complex(1.0, 0.0)}));
    EXPECT_NEAR(std::real(result), -1.0, 1e-8);
    EXPECT_NEAR(std::imag(result), 1.0, 1e-8);
}
