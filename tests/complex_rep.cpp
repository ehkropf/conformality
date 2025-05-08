/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
 *
 * This file is part of Conformality.
 * Conformality is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Conformality is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Conformality. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gtest/gtest.h>
#include "../src/complex.hpp"

TEST(ComplexTest, Construction) {
    Complex z1(3.0, 4.0);
    EXPECT_DOUBLE_EQ(3.0, z1.real());
    EXPECT_DOUBLE_EQ(4.0, z1.imag());
}

TEST(ComplexTest, BasicOperations) {
    Complex z1(3.0, 4.0);
    Complex z2(1.0, 2.0);

    // Addition
    Complex sum = z1 + z2;
    EXPECT_DOUBLE_EQ(4.0, sum.real());
    EXPECT_DOUBLE_EQ(6.0, sum.imag());

    // Multiplication
    Complex product = z1 * z2;
    EXPECT_DOUBLE_EQ(-5.0, product.real());
    EXPECT_DOUBLE_EQ(10.0, product.imag());

    // Division
    Complex quotient = z1 / z2;
    EXPECT_NEAR(2.2, quotient.real(), 1e-10);
    EXPECT_NEAR(-0.4, quotient.imag(), 1e-10);
}

TEST(ComplexTest, PolarForm) {
    Complex z(3.0, 4.0);
    EXPECT_DOUBLE_EQ(5.0, z.abs());
    EXPECT_NEAR(0.9272952180016122, z.arg(), 1e-10);

    Complex fromPolar = ComplexDouble::fromPolar(5.0, 0.9272952180016122);
    EXPECT_NEAR(3.0, fromPolar.real(), 1e-10);
    EXPECT_NEAR(4.0, fromPolar.imag(), 1e-10);
}
