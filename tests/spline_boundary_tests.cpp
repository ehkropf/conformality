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

#include <gtest/gtest.h>
#include "../src/domains/SplineBoundaryComponent.h"
#include "../src/domains/Boundary.h"

#include <cmath>

namespace
{

// Thesis example 1 control points (11 unique + closure)
const std::vector<double> THESIS1_XPTS = {
    1.956140, 1.570175, 0.710526, 0.008772, -0.412281, -1.289474,
    -1.798246, -2.026316, -1.149123, 0.692982, 1.728070, 1.956140
};
const std::vector<double> THESIS1_YPTS = {
    0.043860, 0.500000, 0.657895, 0.815789, 1.429825, 1.605263,
    0.710526, -0.622807, -1.675439, -1.763158, -1.061404, 0.043860
};

// Generate circle control points
std::pair<std::vector<double>, std::vector<double>> circleControlPoints(int n, double radius = 1.0)
{
    std::vector<double> x(n + 1), y(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        double t = 2.0 * M_PI * i / n;
        x[i] = radius * std::cos(t);
        y[i] = radius * std::sin(t);
    }
    return {x, y};
}

} // anonymous namespace

// --- Construction ---

TEST(SplineBoundaryComponent, ConstructsFromThesis1Points)
{
    EXPECT_NO_THROW(SplineBoundaryComponent(THESIS1_XPTS, THESIS1_YPTS));
}

TEST(SplineBoundaryComponent, ConstructsWithRefinement)
{
    EXPECT_NO_THROW(SplineBoundaryComponent(THESIS1_XPTS, THESIS1_YPTS, 256));
}

TEST(SplineBoundaryComponent, ThrowsOnMismatchedSizes)
{
    std::vector<double> x = {0, 1, 0};
    std::vector<double> y = {1, 0};
    EXPECT_THROW(SplineBoundaryComponent(x, y), std::invalid_argument);
}

TEST(SplineBoundaryComponent, ThrowsOnTooFewPoints)
{
    EXPECT_THROW(SplineBoundaryComponent({0, 0}, {1, 1}), std::invalid_argument);
}

// --- totalLength ---

TEST(SplineBoundaryComponent, Thesis1RawTotalLength)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    double tl = spline.totalLength();
    // Raw chordal arc-length of 11 control points is ~11.8
    EXPECT_GT(tl, 10.0);
    EXPECT_LT(tl, 15.0);
}

TEST(SplineBoundaryComponent, Thesis1RefinedTotalLengthGrowsWithRefinement)
{
    SplineBoundaryComponent raw(THESIS1_XPTS, THESIS1_YPTS);
    SplineBoundaryComponent refined(THESIS1_XPTS, THESIS1_YPTS, 256);
    // Refined spline has more control points, so chordal arc-length better
    // approximates true arc-length (should be >= raw chordal length)
    EXPECT_GE(refined.totalLength(), raw.totalLength());
    EXPECT_GT(refined.totalLength(), 10.0);
}

TEST(SplineBoundaryComponent, CircleTotalLengthApprox2Pi)
{
    auto [cx, cy] = circleControlPoints(32);
    SplineBoundaryComponent spline(cx, cy);
    // Chordal arc-length of 32-gon inscribed in unit circle
    // sum of 32 chords of length 2*sin(pi/32) ≈ 6.2735
    double expected_chordal = 32.0 * 2.0 * std::sin(M_PI / 32.0);
    EXPECT_NEAR(spline.totalLength(), expected_chordal, 0.001);
}

// --- Periodicity ---

TEST(SplineBoundaryComponent, EvaluateIsPeriodic)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    double tl = spline.totalLength();
    Complex z0 = spline.evaluate(0.0);
    Complex z_tl = spline.evaluate(tl);
    EXPECT_NEAR(std::abs(z0 - z_tl), 0.0, 1e-10);
}

TEST(SplineBoundaryComponent, EvaluateWrapsNegative)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    double tl = spline.totalLength();
    Complex z_pos = spline.evaluate(1.5);
    Complex z_neg = spline.evaluate(1.5 - tl);
    EXPECT_NEAR(std::abs(z_pos - z_neg), 0.0, 1e-10);
}

// --- Evaluate at control point ---

TEST(SplineBoundaryComponent, EvaluateAtZeroReturnsFirstPoint)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    Complex z = spline.evaluate(0.0);
    EXPECT_NEAR(z.real(), THESIS1_XPTS[0], 1e-10);
    EXPECT_NEAR(z.imag(), THESIS1_YPTS[0], 1e-10);
}

// --- Derivative consistency via finite difference ---

TEST(SplineBoundaryComponent, DerivativeConsistentWithFiniteDifference)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    double tl = spline.totalLength();

    // Test at several points
    for (double frac : {0.1, 0.25, 0.5, 0.75, 0.9})
    {
        double s = frac * tl;
        Complex deriv = spline.evaluateDerivative(s);

        double h = 1e-6;
        Complex fd = (spline.evaluate(s + h) - spline.evaluate(s - h)) / (2.0 * h);

        EXPECT_NEAR(deriv.real(), fd.real(), 1e-4)
            << "Real derivative mismatch at s=" << s;
        EXPECT_NEAR(deriv.imag(), fd.imag(), 1e-4)
            << "Imag derivative mismatch at s=" << s;
    }
}

// --- findParameterization round-trip ---

TEST(SplineBoundaryComponent, FindParameterizationRoundTrip)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    double tl = spline.totalLength();

    for (double frac : {0.0, 0.2, 0.5, 0.8})
    {
        double s_orig = frac * tl;
        Complex z = spline.evaluate(s_orig);
        double s_found = spline.findParameterization(z);
        Complex z_found = spline.evaluate(s_found);

        EXPECT_NEAR(std::abs(z - z_found), 0.0, 1e-8)
            << "Round-trip failed at s=" << s_orig;
    }
}

// --- sample ---

TEST(SplineBoundaryComponent, SampleReturnsCorrectCount)
{
    SplineBoundaryComponent spline(THESIS1_XPTS, THESIS1_YPTS);
    auto samples = spline.sample(256);
    EXPECT_EQ(samples.size(), 256u);
}

// --- Boundary integration ---

TEST(SplineBoundaryComponent, WorksWithBoundaryClass)
{
    auto component = std::make_shared<SplineBoundaryComponent>(THESIS1_XPTS, THESIS1_YPTS, 256);
    Boundary boundary(component);
    EXPECT_EQ(boundary.getNumComponents(), 1u);
    EXPECT_NEAR(boundary.totalLength(), component->totalLength(), 1e-12);
}

// --- Circle with refinement ---

TEST(SplineBoundaryComponent, RefinedCircleApproachesAnalytic)
{
    auto [cx, cy] = circleControlPoints(8);
    // Refine with N=256
    SplineBoundaryComponent spline(cx, cy, 256);
    double tl = spline.totalLength();

    // The refined spline should be closer to the actual unit circle
    double max_error = 0.0;
    for (int j = 0; j < 100; ++j)
    {
        double s = tl * j / 100.0;
        Complex z = spline.evaluate(s);
        double err = std::abs(std::abs(z) - 1.0);
        max_error = std::max(max_error, err);
    }
    // With refinement from 8 control points, should be quite close to circle
    EXPECT_LT(max_error, 0.01);
}
