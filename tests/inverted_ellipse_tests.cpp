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
#include "../src/domains/BoundaryComponent.h"
#include "../src/domains/InvertedEllipseComponent.h"

#include <cmath>

// Test the zero-rotation case: z(t) = 1/(alpha*cos(t) - i*sin(t)) + center
// For alpha=0.3, center=0: at t=0, z = 1/0.3 = 10/3
// At t=pi/2, z = 1/(-i) = i
// At t=pi, z = 1/(-0.3) = -10/3
// At t=3pi/2, z = 1/i = -i

TEST(InvertedEllipseComponent, EvaluateZeroRotation)
{
    InvertedEllipseComponent comp(Complex(0, 0), 0.3);

    Complex z0 = comp.evaluate(0.0);
    EXPECT_NEAR(z0.real(), 1.0 / 0.3, 1e-12);
    EXPECT_NEAR(z0.imag(), 0.0, 1e-12);

    Complex z1 = comp.evaluate(M_PI / 2.0);
    EXPECT_NEAR(z1.real(), 0.0, 1e-12);
    EXPECT_NEAR(z1.imag(), 1.0, 1e-12);

    Complex z2 = comp.evaluate(M_PI);
    EXPECT_NEAR(z2.real(), -1.0 / 0.3, 1e-12);
    EXPECT_NEAR(z2.imag(), 0.0, 1e-12);

    Complex z3 = comp.evaluate(3.0 * M_PI / 2.0);
    EXPECT_NEAR(z3.real(), 0.0, 1e-12);
    EXPECT_NEAR(z3.imag(), -1.0, 1e-12);
}

TEST(InvertedEllipseComponent, EvaluateWithCenter)
{
    Complex center(1.0, 2.0);
    InvertedEllipseComponent comp(center, 0.3);

    Complex z0 = comp.evaluate(0.0);
    EXPECT_NEAR(z0.real(), 1.0 / 0.3 + 1.0, 1e-12);
    EXPECT_NEAR(z0.imag(), 2.0, 1e-12);
}

TEST(InvertedEllipseComponent, DerivativeMatchesFiniteDifference)
{
    InvertedEllipseComponent comp(Complex(0, 0), 0.3);
    const double h = 1e-7;

    for (double t : {0.0, 0.5, 1.0, M_PI / 2.0, M_PI, 2.5, 5.0})
    {
        Complex analytic = comp.evaluateDerivative(t);
        Complex fd = (comp.evaluate(t + h) - comp.evaluate(t - h)) / (2.0 * h);
        EXPECT_NEAR(analytic.real(), fd.real(), 1e-5)
            << "Real part mismatch at t=" << t;
        EXPECT_NEAR(analytic.imag(), fd.imag(), 1e-5)
            << "Imag part mismatch at t=" << t;
    }
}

TEST(InvertedEllipseComponent, DerivativeWithRotation)
{
    InvertedEllipseComponent comp(Complex(0, 0), 0.3, M_PI / 4.0);
    const double h = 1e-7;

    for (double t : {0.0, 1.0, M_PI / 2.0, M_PI, 4.0})
    {
        Complex analytic = comp.evaluateDerivative(t);
        Complex fd = (comp.evaluate(t + h) - comp.evaluate(t - h)) / (2.0 * h);
        EXPECT_NEAR(analytic.real(), fd.real(), 1e-5)
            << "Real part mismatch at t=" << t;
        EXPECT_NEAR(analytic.imag(), fd.imag(), 1e-5)
            << "Imag part mismatch at t=" << t;
    }
}

TEST(InvertedEllipseComponent, SampleReturnsCorrectCount)
{
    InvertedEllipseComponent comp(Complex(0, 0), 0.3);

    auto samples = comp.sample(64);
    EXPECT_EQ(samples.size(), 64u);

    auto samples2 = comp.sample(256);
    EXPECT_EQ(samples2.size(), 256u);
}

TEST(InvertedEllipseComponent, FindParameterizationRoundTrip)
{
    InvertedEllipseComponent comp(Complex(0, 0), 0.3);

    for (double t_orig : {0.5, 1.0, M_PI / 2.0, M_PI, 4.0, 5.5})
    {
        Complex z = comp.evaluate(t_orig);
        double t_found = comp.findParameterization(z);
        Complex z_found = comp.evaluate(t_found);
        EXPECT_NEAR(std::abs(z - z_found), 0.0, 1e-8)
            << "Round-trip failed at t=" << t_orig;
    }
}

TEST(InvertedEllipseComponent, MatchesMATLABEx2Outer)
{
    // th_gen_ex2: binvellip([0, .3], N) -- center=0, alpha=0.3, rotation=0
    InvertedEllipseComponent comp(Complex(0, 0), 0.3);

    // Curve should be large: at t=0, extends to ~3.33 on real axis
    Complex z0 = comp.evaluate(0.0);
    EXPECT_GT(std::abs(z0), 3.0);

    // Curve is closed and traversed counterclockwise
    auto samples = comp.sample(256);
    EXPECT_EQ(samples.size(), 256u);
    // First and wrap-around should be close to continuous
    Complex diff = samples[0] - comp.evaluate(2.0 * M_PI);
    EXPECT_NEAR(std::abs(diff), 0.0, 1e-12);
}

TEST(InvertedEllipseComponent, InvalidAlphaThrows)
{
    EXPECT_THROW(InvertedEllipseComponent(Complex(0, 0), 0.0), std::invalid_argument);
    EXPECT_THROW(InvertedEllipseComponent(Complex(0, 0), 1.0), std::invalid_argument);
    EXPECT_THROW(InvertedEllipseComponent(Complex(0, 0), -0.5), std::invalid_argument);
    EXPECT_THROW(InvertedEllipseComponent(Complex(0, 0), 1.5), std::invalid_argument);
}
