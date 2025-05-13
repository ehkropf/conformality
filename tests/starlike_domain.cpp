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
#include "../src/Domain.hpp"

TEST(StarlikeDomainTest, Construction) {
    Complex center(0.0, 0.0);
    // FIXME This is clearly the wrong radius function.
    auto radiusFunc = [](double _) -> double { return 1.0; }; // Constant radius = circle

    StarlikeDomain domain(center, radiusFunc, false); // Internal domain

    EXPECT_FALSE(domain.isExternalDomain());
    EXPECT_EQ(1, domain.getConnectivity());
    EXPECT_EQ(center, domain.getCenter());
}

TEST(StarlikeDomainTest, ContainsPoint) {
    Complex center(0.0, 0.0);
    // FIXME This is clearly the wrong radius function.
    auto radiusFunc = [](double _) -> double { return 1.0; }; // Circle of radius 1

    StarlikeDomain internalDomain(center, radiusFunc, false);

    // Internal points
    EXPECT_TRUE(internalDomain.contains(Complex(0.0, 0.0)));
    EXPECT_TRUE(internalDomain.contains(Complex(0.5, 0.0)));

    // Boundary point (approximately)
    EXPECT_TRUE(internalDomain.contains(Complex(0.999, 0.0)));

    // External point
    EXPECT_FALSE(internalDomain.contains(Complex(1.5, 0.0)));

    // Now test external domain
    StarlikeDomain externalDomain(center, radiusFunc, true);

    // Internal points (now external to the domain)
    EXPECT_FALSE(externalDomain.contains(Complex(0.0, 0.0)));
    EXPECT_FALSE(externalDomain.contains(Complex(0.5, 0.0)));

    // External points (now inside the domain)
    EXPECT_TRUE(externalDomain.contains(Complex(1.5, 0.0)));
}

TEST(StarlikeDomainTest, RadiusFunction) {
    Complex center(0.0, 0.0);
    // Elliptical-like radius function
    double a = 2.0, b = 1.0;
    auto radiusFunc = [a, b](double angle) -> double {
        return a * b / std::sqrt(b * b * std::cos(angle) * std::cos(angle) +
                               a * a * std::sin(angle) * std::sin(angle));
    };

    StarlikeDomain domain(center, radiusFunc, false);

    // Check radii at different angles
    EXPECT_NEAR(a, domain.getRadius(0.0), 1e-10);           // Along x-axis
    EXPECT_NEAR(b, domain.getRadius(M_PI / 2), 1e-10);      // Along y-axis
}
