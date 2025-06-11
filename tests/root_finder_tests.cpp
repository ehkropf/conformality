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
#include <cmath>
#include <complex>
#include "../src/RootFinder.h"

TEST(RootFinderTest, TernarySearchQuadratic)
{
    // Test with a simple quadratic function: f(x) = (x - 2)^2
    // Minimum should be at x = 2
    auto quadratic = [](double x) -> double
    {
        return (x - 2.0) * (x - 2.0);
    };

    double result = RootFinder::ternarySearch(quadratic, 0.0, 4.0);
    EXPECT_NEAR(2.0, result, 1e-8);
}

TEST(RootFinderTest, TernarySearchCosine)
{
    // Test with cosine function: f(x) = cos(x)
    // Minimum in [0, π] should be at x = π
    auto cosine = [](double x) -> double
    {
        return std::cos(x);
    };

    double result = RootFinder::ternarySearch(cosine, 0.0, M_PI);
    EXPECT_NEAR(M_PI, result, 1e-8);
}

TEST(RootFinderTest, TernarySearchCustomTolerance)
{
    // Test with custom tolerance
    auto quadratic = [](double x) -> double
    {
        return (x - 1.5) * (x - 1.5);
    };

    double result = RootFinder::ternarySearch(quadratic, 0.0, 3.0, 1e-12);
    EXPECT_NEAR(1.5, result, 1e-11);
}

TEST(RootFinderTest, NewtonMethodRealQuadratic)
{
    // Test Newton's method with f(x) = x^2 - 4, f'(x) = 2x
    // Root should be at x = 2 (starting from x = 3)
    auto function = [](double x) -> double
    {
        return x * x - 4.0;
    };
    auto derivative = [](double x) -> double
    {
        return 2.0 * x;
    };

    double result = RootFinder::newton<double>(function, derivative, 3.0);
    EXPECT_NEAR(2.0, result, 1e-8);
}

TEST(RootFinderTest, NewtonMethodRealCubic)
{
    // Test Newton's method with f(x) = x^3 - 2x - 5, f'(x) = 3x^2 - 2
    // Has a root around x ≈ 2.094551
    auto function = [](double x) -> double
    {
        return x * x * x - 2.0 * x - 5.0;
    };
    auto derivative = [](double x) -> double
    {
        return 3.0 * x * x - 2.0;
    };

    double result = RootFinder::newton<double>(function, derivative, 2.0);
    EXPECT_NEAR(2.094551481542327, result, 1e-8);
}

TEST(RootFinderTest, NewtonMethodComplexPolynomial)
{
    // Test Newton's method with complex polynomial f(z) = z^2 + 1
    // f'(z) = 2z, roots at z = ±i
    auto function = [](std::complex<double> z) -> std::complex<double>
    {
        return z * z + 1.0;
    };
    auto derivative = [](std::complex<double> z) -> std::complex<double>
    {
        return 2.0 * z;
    };

    // Starting near +i, should converge to +i
    std::complex<double> initial{0.1, 1.1};
    std::complex<double> result = RootFinder::newton<std::complex<double>>(function, derivative, initial);

    EXPECT_NEAR(0.0, result.real(), 1e-8);
    EXPECT_NEAR(1.0, result.imag(), 1e-8);
}

TEST(RootFinderTest, NewtonMethodComplexCubic)
{
    // Test with f(z) = z^3 - 1, f'(z) = 3z^2
    // Has roots at the cube roots of unity
    auto function = [](std::complex<double> z) -> std::complex<double>
    {
        return z * z * z - 1.0;
    };
    auto derivative = [](std::complex<double> z) -> std::complex<double>
    {
        return 3.0 * z * z;
    };

    // Starting near z = 1, should converge to z = 1
    std::complex<double> initial{1.1, 0.1};
    std::complex<double> result = RootFinder::newton<std::complex<double>>(function, derivative, initial);

    EXPECT_NEAR(1.0, result.real(), 1e-8);
    EXPECT_NEAR(0.0, result.imag(), 1e-8);
}

TEST(RootFinderTest, NewtonMethodCustomTolerance)
{
    // Test Newton's method with custom tolerance
    auto function = [](double x) -> double
    {
        return x * x - 9.0;
    };
    auto derivative = [](double x) -> double
    {
        return 2.0 * x;
    };

    double result = RootFinder::newton<double>(function, derivative, 4.0, 1e-12);
    EXPECT_NEAR(3.0, result, 1e-11);
}

TEST(RootFinderTest, NewtonMethodMaxIterations)
{
    // Test that Newton's method respects max iterations
    // Using a function with a problematic derivative near the root
    auto function = [](double x) -> double
    {
        return x * x * x - 2.0 * x * x + x - 1.0;
    };
    auto derivative = [](double x) -> double
    {
        return 3.0 * x * x - 4.0 * x + 1.0;
    };

    // Should throw convergence error when max iterations reached
    EXPECT_THROW({
        RootFinder::newton<double>(function, derivative, 0.5, 1e-15, 5);
    }, RootFinder::ConvergenceError);
}

TEST(RootFinderTest, TernarySearchEdgeCases)
{
    // Test when minimum is at the boundary
    auto function = [](double x) -> double
    {
        return x; // linear function, minimum at left boundary
    };

    double result = RootFinder::ternarySearch(function, 1.0, 3.0);
    EXPECT_NEAR(1.0, result, 1e-8);
}

TEST(RootFinderTest, NewtonMethodConvergenceToZero)
{
    // Test finding a root where the function value is exactly zero
    auto function = [](double x) -> double
    {
        return (x - 1.0) * (x - 1.0) * (x - 1.0); // Triple root at x = 1
    };
    auto derivative = [](double x) -> double
    {
        return 3.0 * (x - 1.0) * (x - 1.0);
    };

    double result = RootFinder::newton<double>(function, derivative, 1.5);
    EXPECT_NEAR(1.0, result, 1e-3); // Slower convergence for multiple roots
}
