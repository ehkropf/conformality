/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
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
#include "../src/Complex.h"

TEST(ComplexTest, Construction)
{
    Complex z1(3.0, 4.0);
    EXPECT_DOUBLE_EQ(3.0, z1.real());
    EXPECT_DOUBLE_EQ(4.0, z1.imag());
}

TEST(ComplexTest, DefaultConstructor)
{
    Complex z;
    EXPECT_DOUBLE_EQ(0.0, z.real());
    EXPECT_DOUBLE_EQ(0.0, z.imag());
}

TEST(ComplexTest, ConstructFromComplex)
{
    std::complex<double> stdz(2.0, 3.0);
    Complex z(stdz);
    EXPECT_DOUBLE_EQ(2.0, z.real());
    EXPECT_DOUBLE_EQ(3.0, z.imag());
}

TEST(ComplexTest, BasicOperations)
{
    Complex z1(3.0, 4.0);
    Complex z2(1.0, 2.0);

    // Addition
    Complex sum = z1 + z2;
    EXPECT_DOUBLE_EQ(4.0, sum.real());
    EXPECT_DOUBLE_EQ(6.0, sum.imag());

    // Complex
    Complex diff = z1 - z2;
    EXPECT_DOUBLE_EQ(2.0, diff.real());
    EXPECT_DOUBLE_EQ(2.0, diff.real());

    // Multiplication
    Complex product = z1 * z2;
    EXPECT_DOUBLE_EQ(-5.0, product.real());
    EXPECT_DOUBLE_EQ(10.0, product.imag());

    // Division
    Complex quotient = z1 / z2;
    EXPECT_NEAR(2.2, quotient.real(), 1e-10);
    EXPECT_NEAR(-0.4, quotient.imag(), 1e-10);
}

TEST(ComplexTest, GetValue)
{
    Complex z(1.0, 2.0);
    std::complex<double> stdz = z.getValue();
    EXPECT_DOUBLE_EQ(1.0, std::real(stdz));
    EXPECT_DOUBLE_EQ(2.0, std::imag(stdz));

    // Test non-const getValue
    Complex z2(3.0, 4.0);
    z2.getValue() = std::complex<double>(5.0, 6.0);
    EXPECT_DOUBLE_EQ(5.0, z2.real());
    EXPECT_DOUBLE_EQ(6.0, z2.imag());
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
    EXPECT_DOUBLE_EQ(5.0, z.abs());
    EXPECT_NEAR(0.9272952180016122, z.arg(), 1e-10);

    Complex fromPolar = ComplexDouble::fromPolar(5.0, 0.9272952180016122);
    EXPECT_NEAR(3.0, fromPolar.real(), 1e-10);
    EXPECT_NEAR(4.0, fromPolar.imag(), 1e-10);
}

TEST(ComplexTest, DifferentTemplateType)
{
    // Test with float instead of double
    Complex<std::complex<float>> zf(1.0f, 2.0f);
    EXPECT_FLOAT_EQ(1.0f, zf.real());
    EXPECT_FLOAT_EQ(2.0f, zf.imag());

    // Basic operations
    Complex<std::complex<float>> zf2(2.0f, 3.0f);
    Complex<std::complex<float>> sum = zf + zf2;
    EXPECT_FLOAT_EQ(3.0f, sum.real());
    EXPECT_FLOAT_EQ(5.0f, sum.imag());
}
