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
#include "../src/Types.h"

TEST(ComplexTest, Construction)
{
    Complex z1(3.0, 4.0);
    EXPECT_DOUBLE_EQ(3.0, std::real(z1));
    EXPECT_DOUBLE_EQ(4.0, std::imag(z1));
}

TEST(ComplexTest, DefaultConstructor)
{
    Complex z;
    EXPECT_DOUBLE_EQ(0.0, std::real(z));
    EXPECT_DOUBLE_EQ(0.0, std::imag(z));
}

TEST(ComplexTest, ConstructFromComplex)
{
    std::complex<double> stdz(2.0, 3.0);
    Complex z = stdz;
    EXPECT_DOUBLE_EQ(2.0, std::real(z));
    EXPECT_DOUBLE_EQ(3.0, std::imag(z));
}

TEST(ComplexTest, BasicOperations)
{
    Complex z1(3.0, 4.0);
    Complex z2(1.0, 2.0);

    // Addition
    Complex sum = z1 + z2;
    EXPECT_DOUBLE_EQ(4.0, std::real(sum));
    EXPECT_DOUBLE_EQ(6.0, std::imag(sum));

    // Subtraction
    Complex diff = z1 - z2;
    EXPECT_DOUBLE_EQ(2.0, std::real(diff));
    EXPECT_DOUBLE_EQ(2.0, std::imag(diff));

    // Multiplication
    Complex product = z1 * z2;
    EXPECT_DOUBLE_EQ(-5.0, std::real(product));
    EXPECT_DOUBLE_EQ(10.0, std::imag(product));

    // Division
    Complex quotient = z1 / z2;
    EXPECT_NEAR(2.2, std::real(quotient), 1e-10);
    EXPECT_NEAR(-0.4, std::imag(quotient), 1e-10);
}

TEST(ComplexTest, Assignment)
{
    Complex z(1.0, 2.0);
    std::complex<double> stdz = z;
    EXPECT_DOUBLE_EQ(1.0, std::real(stdz));
    EXPECT_DOUBLE_EQ(2.0, std::imag(stdz));

    // Test assignment
    Complex z2(3.0, 4.0);
    z2 = std::complex<double>(5.0, 6.0);
    EXPECT_DOUBLE_EQ(5.0, std::real(z2));
    EXPECT_DOUBLE_EQ(6.0, std::imag(z2));
}

TEST(ComplexTest, EqualityOperator)
{
    Complex z1(1.0, 2.0);
    Complex z2(1.0, 2.0);
    Complex z3(3.0, 4.0);

    EXPECT_TRUE(z1 == z2);
    EXPECT_FALSE(z1 == z3);
}

TEST(ComplexTest, PolarForm)
{
    Complex z(3.0, 4.0);
    EXPECT_DOUBLE_EQ(5.0, std::abs(z));
    EXPECT_NEAR(0.9272952180016122, std::arg(z), 1e-10);

    Complex fromPolar = std::polar(5.0, 0.9272952180016122);
    EXPECT_NEAR(3.0, std::real(fromPolar), 1e-10);
    EXPECT_NEAR(4.0, std::imag(fromPolar), 1e-10);
}

TEST(ComplexTest, FloatComplex)
{
    // Test with float complex numbers
    std::complex<float> zf(1.0f, 2.0f);
    EXPECT_FLOAT_EQ(1.0f, std::real(zf));
    EXPECT_FLOAT_EQ(2.0f, std::imag(zf));

    // Basic operations
    std::complex<float> zf2(2.0f, 3.0f);
    std::complex<float> sum = zf + zf2;
    EXPECT_FLOAT_EQ(3.0f, std::real(sum));
    EXPECT_FLOAT_EQ(5.0f, std::imag(sum));
}
