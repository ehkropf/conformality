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

#include "../src/methods/MCSCUnboundedObjectiveFunction.h"

#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace
{

// example_driver_ext.m's 3-polygon unbounded example: three triangles as the target, with a
// matching 3-circle initial guess (equilateral-ish circles at the prevertex angles given in the
// MATLAB source).
MCSCPolygonalDomain makeExampleDriverTargetPolygon()
{
    std::vector<Complex> triangle0 = {Complex(0.0, 0.0), Complex(-1.5, -2.0), Complex(1.5, -2.0)};
    std::vector<Complex> triangle1 = {
        Complex(-1.0, 1.0), Complex(-1.5, 1.0 + std::sqrt(3.0) / 2.0), Complex(-2.0, 1.0)};
    std::vector<Complex> triangle2 = {
        Complex(1.0, 1.0), Complex(2.0, 1.0), Complex(1.5, 1.0 + std::sqrt(3.0) / 2.0)};
    return MCSCPolygonalDomain({triangle0, triangle1, triangle2}, /*isUnboundedDomain=*/true);
}

MCSCCircleDomain makeExampleDriverInitialGuess()
{
    const double twoPiThird = 2.0 * M_PI / 3.0;
    // example_driver_ext.m's angles for circle 1 ([3*pi/2; pi/4; 3*pi/4]) are given in
    // MATLAB's own unconstrained-parameter convention (only the first angle need be in
    // [0, 2*pi); the rest are implicitly unwound to stay strictly increasing) -- unwind them
    // here since MCSCCircleDomain requires an explicitly strictly-increasing input sequence.
    return MCSCCircleDomain(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, twoPiThird, 2.0 * twoPiThird}},
        {Complex(2.4, 1.3), 0.37, {3.0 * M_PI / 2.0, M_PI / 4.0 + 2.0 * M_PI, 3.0 * M_PI / 4.0 + 2.0 * M_PI}},
        {Complex(2.4, -1.3), 0.37, {M_PI / 2.0, 5.0 * M_PI / 4.0, 7.0 * M_PI / 4.0}},
    });
}

} // namespace

TEST(MCSCUnboundedObjectiveFunctionTest, EvaluateReturnsFiniteVectorOfExpectedLength)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto circle = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, circle);

    Eigen::VectorXd Xu = circle.toUnconstrained();
    Eigen::VectorXd F = objective.evaluate(Xu);

    EXPECT_EQ(F.size(), Xu.size());
    EXPECT_TRUE(F.allFinite());
}

TEST(MCSCUnboundedObjectiveFunctionTest, ConstructorPropagatesConnectivityMismatch)
{
    auto circle = makeExampleDriverInitialGuess();  // connectivity 3
    MCSCPolygonalDomain twoComponentPolygon(
        {{Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.5, 1.0)},
         {Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(5.5, 1.0)}},
        /*isUnboundedDomain=*/true);

    EXPECT_THROW(MCSCUnboundedObjectiveFunction(twoComponentPolygon, circle), std::invalid_argument);
}

TEST(MCSCUnboundedObjectiveFunctionTest, ConstructorPropagatesVertexCountMismatch)
{
    auto circle = makeExampleDriverInitialGuess();  // 3 prevertices per circle
    std::vector<Complex> squareA = {
        Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 1.0), Complex(0.0, 1.0)};
    std::vector<Complex> triangleB = {Complex(5.0, 0.0), Complex(6.0, 0.0), Complex(5.5, 1.0)};
    std::vector<Complex> triangleC = {Complex(-5.0, 0.0), Complex(-6.0, 0.0), Complex(-5.5, 1.0)};
    MCSCPolygonalDomain mismatched({squareA, triangleB, triangleC}, /*isUnboundedDomain=*/true);

    EXPECT_THROW(MCSCUnboundedObjectiveFunction(mismatched, circle), std::invalid_argument);
}

TEST(MCSCUnboundedObjectiveFunctionTest, EvaluateThrowsOnWrongLengthXu)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto circle = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, circle);

    Eigen::VectorXd tooShort(1);
    tooShort << 0.0;

    EXPECT_THROW(objective.evaluate(tooShort), std::invalid_argument);
}

TEST(MCSCUnboundedObjectiveFunctionTest, EvaluateThrowsOnNonFiniteXu)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto circle = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, circle);

    Eigen::VectorXd Xu = circle.toUnconstrained();
    Xu[0] = std::numeric_limits<double>::quiet_NaN();

    EXPECT_THROW(objective.evaluate(Xu), std::invalid_argument);
}

TEST(MCSCUnboundedObjectiveFunctionTest, SingleCircleConnectivityIsRejectedByReflectionMachinery)
{
    // m == 1 is not a supported configuration: the reflection method (mcsc::
    // reflectCircleSequence, landed in #163) needs at least one other circle to reflect through,
    // and throws accordingly. Not a defect introduced by this PR -- documenting the boundary
    // behavior here rather than asserting m==1 "works", since it structurally cannot.
    std::vector<Complex> triangle = {Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.5, 1.0)};
    MCSCPolygonalDomain polygon({triangle}, /*isUnboundedDomain=*/true);
    MCSCCircleDomain circle(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0}},
    });

    EXPECT_THROW(MCSCUnboundedObjectiveFunction(polygon, circle), std::invalid_argument);
}

TEST(MCSCUnboundedObjectiveFunctionTest, EvaluateIsDeterministic)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto circle = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, circle);

    Eigen::VectorXd Xu = circle.toUnconstrained();
    Eigen::VectorXd F1 = objective.evaluate(Xu);
    Eigen::VectorXd F2 = objective.evaluate(Xu);

    for (int i = 0; i < F1.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(F1[i], F2[i]);
    }
}

TEST(MCSCUnboundedObjectiveFunctionTest, EvaluateUpdatesWorkingCircleDomain)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto circle = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, circle);

    Eigen::VectorXd Xu = circle.toUnconstrained();
    Xu[0] += 0.05;  // perturb log(r_1)
    objective.evaluate(Xu);

    EXPECT_NEAR(objective.workingCircleDomain().getRadius(1), circle.getRadius(1) * std::exp(0.05), 1e-12);
}

// End-to-end regression test: solve the example_driver_ext.m parameter problem via continuation
// and verify convergence + geometric sanity, without needing MATLAB-generated reference numbers
// (full numeric comparison against MATLAB is #169's job).
TEST(MCSCUnboundedObjectiveFunctionTest, SolvesExampleDriverExtParameterProblem)
{
    auto polygon = makeExampleDriverTargetPolygon();
    auto initialGuess = makeExampleDriverInitialGuess();
    MCSCUnboundedObjectiveFunction objective(polygon, initialGuess);

    MCSCCircleDomain solved = objective.solve();

    Eigen::VectorXd Xu = solved.toUnconstrained();
    Eigen::VectorXd F = objective.evaluate(Xu);
    EXPECT_LT(F.lpNorm<Eigen::Infinity>(), 1e-8);

    // Independent sanity check: recompute each polygon's actual side lengths from the solved
    // circle domain (via the integrand + quadrature directly, not the objective's own residual
    // bookkeeping) and compare against the target polygon's true side lengths, up to the same
    // scaling constant A used internally.
    MCSCReflectionIntegrand integrand(polygon, solved, 6);
    GaussJacobiQuadrature::Integrand fprime = [&integrand](const Complex& z)
    { return integrand.evalFPrime(z); };

    std::vector<double> betaValues;
    for (int j = 0; j < polygon.getConnectivity(); ++j)
    {
        for (double a : polygon.getAlpha(j))
        {
            betaValues.push_back(1.0 - a);
        }
    }
    GaussJacobiQuadrature quadrature(betaValues, 12);

    const auto& t0 = solved.getPrevertexAngles(0);
    const double right0 = (t0[0] > t0[1]) ? t0[1] + 2.0 * M_PI : t0[1];
    const Complex Q12 =
        quadrature.integrateArc(t0[0], 0, right0, 1, solved.getCenter(0), solved.getRadius(0), fprime, {});
    const auto& w0 = polygon.getVertices(0);
    const Complex A = (w0[1] - w0[0]) / Q12;

    int offset = 0;
    for (int j = 0; j < polygon.getConnectivity(); ++j)
    {
        const auto& t = solved.getPrevertexAngles(j);
        const auto& w = polygon.getVertices(j);
        const int Kj = static_cast<int>(t.size());
        for (int k = 0; k < Kj; ++k)
        {
            const int kNext = (k + 1) % Kj;
            const double left = t[k];
            const double right = (kNext == 0) ? ((t[0] > left) ? t[0] : t[0] + 2.0 * M_PI)
                                               : ((t[kNext] > left) ? t[kNext] : t[kNext] + 2.0 * M_PI);
            const Complex sideLength = A
                * quadrature.integrateArc(
                    left, offset + k, right, offset + kNext, solved.getCenter(j), solved.getRadius(j), fprime, {});
            const double expected = std::abs(w[kNext] - w[k]);
            EXPECT_NEAR(std::abs(sideLength), expected, 1e-6 * std::max(1.0, expected));
        }
        offset += Kj;
    }
}
