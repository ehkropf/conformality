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
#include "../src/Domain.h"

TEST(EllipticalDomainTest, Construction) {
    double a = 2.0, b = 1.0;
    double rotation = 0.0;
    Complex center(0.0, 0.0);

    EllipticalDomain domain(a, b, rotation, center, false);

    EXPECT_FALSE(domain.isExternalDomain());
    EXPECT_EQ(1, domain.getConnectivity());
    EXPECT_EQ(center, domain.getCenter());
    EXPECT_DOUBLE_EQ(std::sqrt(1.0 - (b*b)/(a*a)), domain.getEccentricity());
}

TEST(EllipticalDomainTest, ContainsPoint) {
    double a = 2.0, b = 1.0;
    double rotation = 0.0;
    Complex center(0.0, 0.0);

    EllipticalDomain domain(a, b, rotation, center, false);

    // Inside points
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.0)));
    EXPECT_TRUE(domain.contains(Complex(1.0, 0.0)));
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.5)));

    // Boundary points (approximately)
    EXPECT_TRUE(domain.contains(Complex(1.99, 0.0)));
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.99)));

    // Outside points
    EXPECT_FALSE(domain.contains(Complex(2.5, 0.0)));
    EXPECT_FALSE(domain.contains(Complex(0.0, 1.5)));

    // Point that would be inside a circle of radius a, but outside the ellipse
    EXPECT_FALSE(domain.contains(Complex(1.5, 0.75)));
}

TEST(EllipticalDomainTest, RotatedEllipse) {
    double a = 2.0, b = 1.0;
    double rotation = M_PI / 4;  // 45 degrees
    Complex center(0.0, 0.0);

    EllipticalDomain domain(a, b, rotation, center, false);

    // Test points along the rotated axes
    double factor = std::sqrt(2.0) / 2.0;  // cos(45°) = sin(45°) = √2/2

    // Points along the semi-major axis (rotated)
    Complex majorAxis(a * factor, a * factor);  // Rotated a units along 45°
    EXPECT_TRUE(domain.contains(majorAxis * 0.99));  // Just inside
    EXPECT_FALSE(domain.contains(majorAxis * 1.01));  // Just outside

    // Points along the semi-minor axis (rotated)
    Complex minorAxis(-b * factor, b * factor);  // Rotated b units along 135°
    EXPECT_TRUE(domain.contains(minorAxis * 0.99));  // Just inside
    EXPECT_FALSE(domain.contains(minorAxis * 1.01));  // Just outside
}
