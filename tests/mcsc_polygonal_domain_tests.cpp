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
#include "../src/domains/MCSCPolygonalDomain.h"

#include <cmath>

namespace
{
    // Unit square, counterclockwise, vertices at the corners.
    std::vector<Complex> unitSquareCcw()
    {
        return {
            Complex(0.0, 0.0),
            Complex(1.0, 0.0),
            Complex(1.0, 1.0),
            Complex(0.0, 1.0),
        };
    }

    // Same square, clockwise orientation.
    std::vector<Complex> unitSquareCw()
    {
        auto v = unitSquareCcw();
        std::reverse(v.begin(), v.end());
        return v;
    }

    // Small square hole, clockwise as seen from outside (domain on the left convention
    // means an interior hole boundary is traversed clockwise for a bounded domain).
    std::vector<Complex> smallSquareCw(Complex center, double halfSide)
    {
        return {
            center + Complex(-halfSide, -halfSide),
            center + Complex(-halfSide, halfSide),
            center + Complex(halfSide, halfSide),
            center + Complex(halfSide, -halfSide),
        };
    }
}

TEST(MCSCPolygonalDomainTest, SquareHasFourRightAngles)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    EXPECT_EQ(1, domain.getConnectivity());
    EXPECT_FALSE(domain.isUnbounded());
    ASSERT_EQ(4, domain.vertexCount(0));

    const auto& alpha = domain.getAlpha(0);
    ASSERT_EQ(4u, alpha.size());
    for (double a : alpha)
    {
        EXPECT_NEAR(0.5, a, 1e-12);
    }
}

TEST(MCSCPolygonalDomainTest, VertexListRoundTrips)
{
    auto square = unitSquareCcw();
    MCSCPolygonalDomain domain({square});

    const auto& vertices = domain.getVertices(0);
    ASSERT_EQ(square.size(), vertices.size());
    for (size_t i = 0; i < square.size(); ++i)
    {
        EXPECT_EQ(square[i], vertices[i]);
    }
}

TEST(MCSCPolygonalDomainTest, ClockwisePolygonIsReoriented)
{
    MCSCPolygonalDomain domain({unitSquareCw()});

    // calc_angles should detect the wrong orientation and flip vertices + angles so the
    // stored data always matches the "domain on the left" convention.
    const auto& vertices = domain.getVertices(0);
    const auto& alpha = domain.getAlpha(0);

    ASSERT_EQ(4u, vertices.size());
    ASSERT_EQ(4u, alpha.size());

    // After the flip, the stored vertex order should match the CCW square traversal
    // reversed relative to the original CW input, i.e. equal to the reverse of the input.
    auto expected = unitSquareCw();
    std::reverse(expected.begin(), expected.end());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_EQ(expected[i], vertices[i]);
    }

    // Note: calc_angles's orientation flip reorders vertices/angles but does not
    // recompute angle values (matching polygon.m) -- a CW-wound square still has raw
    // turning angles of 1.5 (reflex, from the wrong-orientation perspective).
    for (double a : alpha)
    {
        EXPECT_NEAR(1.5, a, 1e-12);
    }
}

TEST(MCSCPolygonalDomainTest, TriangleAnglesSumCorrectly)
{
    std::vector<Complex> triangle = {
        Complex(0.0, 0.0),
        Complex(1.0, 0.0),
        Complex(0.0, 1.0),
    };
    MCSCPolygonalDomain domain({triangle});

    const auto& alpha = domain.getAlpha(0);
    ASSERT_EQ(3u, alpha.size());

    double sum = 0.0;
    for (double a : alpha)
    {
        sum += (1.0 - a);
    }
    EXPECT_NEAR(2.0, sum, 1e-10);
}

TEST(MCSCPolygonalDomainTest, InvalidPolygonThrows)
{
    // A repeated vertex produces a zero-length edge, which makes the turning-angle
    // computation at that vertex degenerate (angle of a zero vector is taken as 0),
    // breaking the sum(1-alpha) == integer invariant that calc_angles checks.
    std::vector<Complex> degenerate = {
        Complex(0.0, 0.0),
        Complex(0.0, 0.0),
        Complex(1.0, 0.0),
        Complex(0.0, 1.0),
    };

    EXPECT_THROW(MCSCPolygonalDomain({degenerate}), std::invalid_argument);
}

TEST(MCSCPolygonalDomainTest, TooFewVerticesThrows)
{
    std::vector<Complex> segment = {Complex(0.0, 0.0), Complex(1.0, 0.0)};
    EXPECT_THROW(MCSCPolygonalDomain({segment}), std::invalid_argument);
}

TEST(MCSCPolygonalDomainTest, EmptyComponentListThrows)
{
    EXPECT_THROW(MCSCPolygonalDomain({}), std::invalid_argument);
}

TEST(MCSCPolygonalDomainTest, MultiplyConnectedConstruction)
{
    auto outer = unitSquareCcw();
    // Scale outer square up so a small hole comfortably fits inside.
    for (auto& v : outer)
    {
        v *= 4.0;
    }
    auto hole = smallSquareCw(Complex(2.0, 2.0), 0.5);

    MCSCPolygonalDomain domain({outer, hole});

    EXPECT_EQ(2, domain.getConnectivity());
    EXPECT_EQ(4, domain.vertexCount(0));
    EXPECT_EQ(4, domain.vertexCount(1));
}

TEST(MCSCPolygonalDomainTest, ContainsSimplyConnectedSquare)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    EXPECT_TRUE(domain.contains(Complex(0.5, 0.5)));
    EXPECT_FALSE(domain.contains(Complex(2.0, 2.0)));
    EXPECT_FALSE(domain.contains(Complex(-0.5, 0.5)));
}

TEST(MCSCPolygonalDomainTest, ContainsUnboundedComplement)
{
    MCSCPolygonalDomain domain({unitSquareCcw()}, /*isUnboundedDomain=*/true);

    EXPECT_FALSE(domain.contains(Complex(0.5, 0.5)));
    EXPECT_TRUE(domain.contains(Complex(2.0, 2.0)));
}

TEST(MCSCPolygonalDomainTest, ContainsRespectsHole)
{
    auto outer = unitSquareCcw();
    for (auto& v : outer)
    {
        v *= 4.0;
    }
    auto hole = smallSquareCw(Complex(2.0, 2.0), 0.5);

    MCSCPolygonalDomain domain({outer, hole});

    EXPECT_TRUE(domain.contains(Complex(0.5, 0.5)));   // Inside outer, outside hole.
    EXPECT_FALSE(domain.contains(Complex(2.0, 2.0)));  // Inside hole.
    EXPECT_FALSE(domain.contains(Complex(10.0, 10.0))); // Outside outer boundary entirely.
}

TEST(MCSCPolygonalDomainTest, TransformBoundaryRecomputesAngles)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    // Scale and translate: should remain a valid square with the same interior angles.
    domain.transformBoundary([](const Complex& z) { return 2.0 * z + Complex(1.0, 1.0); });

    const auto& vertices = domain.getVertices(0);
    ASSERT_EQ(4u, vertices.size());
    EXPECT_EQ(Complex(1.0, 1.0), vertices[0]);
    EXPECT_EQ(Complex(3.0, 1.0), vertices[1]);
    EXPECT_EQ(Complex(3.0, 3.0), vertices[2]);
    EXPECT_EQ(Complex(1.0, 3.0), vertices[3]);

    for (double a : domain.getAlpha(0))
    {
        EXPECT_NEAR(0.5, a, 1e-12);
    }
}

TEST(MCSCPolygonalDomainTest, TransformBoundaryInvalidatingThrows)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    // Collapse only the first vertex onto the second: the resulting zero-length edge
    // makes the turning-angle computation degenerate there, breaking the angle-sum
    // invariant (unlike collapsing all vertices together, which stays self-consistent).
    EXPECT_THROW(
        domain.transformBoundary(
            [](const Complex& z) { return (z == Complex(0.0, 0.0)) ? Complex(1.0, 0.0) : z; }),
        std::invalid_argument);
}

TEST(MCSCPolygonalDomainTest, OutOfRangeComponentIndexThrows)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    EXPECT_THROW(domain.getVertices(1), std::invalid_argument);
    EXPECT_THROW(domain.getAlpha(-1), std::invalid_argument);
    EXPECT_THROW(domain.vertexCount(5), std::invalid_argument);
}
