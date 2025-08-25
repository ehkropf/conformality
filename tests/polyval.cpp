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
#include "../src/numerics/Polyval.h"
#include "../src/core/Types.h"
#include <vector>

class PolyvalTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Setup test data
    }
};

TEST_F(PolyvalTest, SimplePolynomial)
{
    // Test polynomial 2x^2 + 3x + 1 at x = 2
    // Coefficients in descending order: [2, 3, 1]
    // Expected result: 2*4 + 3*2 + 1 = 8 + 6 + 1 = 15
    
    std::vector<double> coeffs{2.0, 3.0, 1.0};
    double x = 2.0;
    double result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, 15.0, 1e-12);
}

TEST_F(PolyvalTest, ConstantPolynomial)
{
    // Test constant polynomial 5 at x = 10
    // Coefficients: [5]
    // Expected result: 5
    
    std::vector<double> coeffs{5.0};
    double x = 10.0;
    double result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, 5.0, 1e-12);
}

TEST_F(PolyvalTest, LinearPolynomial)
{
    // Test linear polynomial 3x + 7 at x = 4
    // Coefficients in descending order: [3, 7]
    // Expected result: 3*4 + 7 = 12 + 7 = 19
    
    std::vector<double> coeffs{3.0, 7.0};
    double x = 4.0;
    double result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, 19.0, 1e-12);
}

TEST_F(PolyvalTest, ComplexPolynomial)
{
    // Test complex polynomial (1+i)z^2 + (2-i)z + 3 at z = 1+i
    // Coefficients: [(1+i), (2-i), 3]
    // z^2 = (1+i)^2 = 1 + 2i - 1 = 2i
    // Expected result: (1+i)*2i + (2-i)*(1+i) + 3
    //                = 2i - 2 + (2+2i-i+1) + 3
    //                = 2i - 2 + 3 + i + 3
    //                = 4 + 3i
    
    std::vector<Complex> coeffs{Complex(1, 1), Complex(2, -1), Complex(3, 0)};
    Complex z{1, 1};
    Complex result = polyval(coeffs, z);
    
    EXPECT_NEAR(result.real(), 4.0, 1e-12);
    EXPECT_NEAR(result.imag(), 3.0, 1e-12);
}

TEST_F(PolyvalTest, ZeroEvaluation)
{
    // Test polynomial 2x^3 + 3x^2 + 4x + 5 at x = 0
    // Expected result: 5 (constant term)
    
    std::vector<double> coeffs{2.0, 3.0, 4.0, 5.0};
    double x = 0.0;
    double result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, 5.0, 1e-12);
}

TEST_F(PolyvalTest, NegativeCoefficients)
{
    // Test polynomial -x^2 + 2x - 1 at x = 3
    // Expected result: -9 + 6 - 1 = -4
    
    std::vector<double> coeffs{-1.0, 2.0, -1.0};
    double x = 3.0;
    double result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, -4.0, 1e-12);
}

// Death test for empty coefficient vector (when assertions are enabled)
#ifdef DEBUG
TEST_F(PolyvalTest, EmptyCoefficientsDeathTest)
{
    std::vector<double> empty_coeffs{};
    double x = 1.0;
    
    EXPECT_DEATH(polyval(empty_coeffs, x), "coefficient vector cannot be empty");
}
#endif

// Test that the concept works with different types
TEST_F(PolyvalTest, ConceptWorksWithFloat)
{
    std::vector<float> coeffs{2.0f, 3.0f, 1.0f};
    float x = 2.0f;
    float result = polyval(coeffs, x);
    
    EXPECT_NEAR(result, 15.0f, 1e-6f);
}

TEST_F(PolyvalTest, ConceptWorksWithInt)
{
    std::vector<int> coeffs{2, 3, 1};
    int x = 2;
    int result = polyval(coeffs, x);
    
    EXPECT_EQ(result, 15);
}