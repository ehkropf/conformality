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
#include "../src/examples/BoundaryHelpers.h"

#include <cmath>

using namespace conformality::examples;

// --- Circle helpers ---

TEST(BoundaryHelpers, CircleEvaluatesCorrectly)
{
    auto boundary = createCircularBoundary(Complex(1.0, 2.0), 0.5);

    Complex z0 = boundary->evaluate(0.0);
    EXPECT_NEAR(z0.real(), 1.5, 1e-12);
    EXPECT_NEAR(z0.imag(), 2.0, 1e-12);

    Complex z1 = boundary->evaluate(M_PI / 2.0);
    EXPECT_NEAR(z1.real(), 1.0, 1e-12);
    EXPECT_NEAR(z1.imag(), 2.5, 1e-12);
}

TEST(BoundaryHelpers, CircleDerivativeCorrect)
{
    auto boundary = createCircularBoundary(Complex(0, 0), 1.0);

    Complex d0 = boundary->evaluateDerivative(0.0);
    EXPECT_NEAR(d0.real(), 0.0, 1e-12);
    EXPECT_NEAR(d0.imag(), 1.0, 1e-12);
}

// --- Ellipse helpers ---

TEST(BoundaryHelpers, EllipseCardinalPoints)
{
    // Axis-aligned ellipse: center=0, a=2, b=1
    auto boundary = createEllipseBoundary(Complex(0, 0), 2.0, 1.0);

    Complex z0 = boundary->evaluate(0.0);
    EXPECT_NEAR(z0.real(), 2.0, 1e-12);
    EXPECT_NEAR(z0.imag(), 0.0, 1e-12);

    Complex z1 = boundary->evaluate(M_PI / 2.0);
    EXPECT_NEAR(z1.real(), 0.0, 1e-12);
    EXPECT_NEAR(z1.imag(), 1.0, 1e-12);

    Complex z2 = boundary->evaluate(M_PI);
    EXPECT_NEAR(z2.real(), -2.0, 1e-12);
    EXPECT_NEAR(z2.imag(), 0.0, 1e-12);
}

TEST(BoundaryHelpers, EllipseWithCenter)
{
    auto boundary = createEllipseBoundary(Complex(1.0, 2.0), 3.0, 1.5);

    Complex z0 = boundary->evaluate(0.0);
    EXPECT_NEAR(z0.real(), 4.0, 1e-12);
    EXPECT_NEAR(z0.imag(), 2.0, 1e-12);
}

TEST(BoundaryHelpers, EllipseWithRotation)
{
    // Ellipse with pi/2 rotation: semi-major along y-axis
    auto boundary = createEllipseBoundary(Complex(0, 0), 2.0, 1.0, M_PI / 2.0);

    Complex z0 = boundary->evaluate(0.0);
    EXPECT_NEAR(z0.real(), 0.0, 1e-12);
    EXPECT_NEAR(z0.imag(), 2.0, 1e-12);
}

TEST(BoundaryHelpers, EllipseDerivativeMatchesFiniteDifference)
{
    auto boundary = createEllipseBoundary(Complex(0, 0), 2.0, 1.0, M_PI / 4.0);
    const double h = 1e-7;

    for (double t : {0.0, 1.0, M_PI / 2.0, M_PI, 4.0})
    {
        Complex analytic = boundary->evaluateDerivative(t);
        Complex fd = (boundary->evaluate(t + h) - boundary->evaluate(t - h)) / (2.0 * h);
        EXPECT_NEAR(analytic.real(), fd.real(), 1e-5) << "at t=" << t;
        EXPECT_NEAR(analytic.imag(), fd.imag(), 1e-5) << "at t=" << t;
    }
}

// --- Inverted ellipse helper ---

TEST(BoundaryHelpers, InvertedEllipseWrapsComponent)
{
    auto boundary = createInvertedEllipseBoundary(Complex(0, 0), 0.3);

    Complex z0 = boundary->evaluate(0.0);
    EXPECT_NEAR(z0.real(), 1.0 / 0.3, 1e-12);
    EXPECT_NEAR(z0.imag(), 0.0, 1e-12);
}

// --- Input validation ---

TEST(BoundaryHelpers, CircleInvalidRadiusThrows)
{
    EXPECT_THROW(createCircularBoundary(Complex(0, 0), 0.0), std::invalid_argument);
    EXPECT_THROW(createCircularBoundary(Complex(0, 0), -1.0), std::invalid_argument);
}

TEST(BoundaryHelpers, EllipseInvalidSemiAxesThrows)
{
    EXPECT_THROW(createEllipseBoundary(Complex(0, 0), 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(createEllipseBoundary(Complex(0, 0), 1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(createEllipseBoundary(Complex(0, 0), -1.0, 1.0), std::invalid_argument);
}

TEST(BoundaryHelpers, InvertedEllipseInvalidAlphaPropagates)
{
    EXPECT_THROW(createInvertedEllipseBoundary(Complex(0, 0), 0.0), std::invalid_argument);
    EXPECT_THROW(createInvertedEllipseBoundary(Complex(0, 0), 1.0), std::invalid_argument);
}
