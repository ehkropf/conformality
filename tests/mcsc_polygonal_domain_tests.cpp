/*
 * Copyright © 2026, Everett Kropf (ehkropf@gmail.com)
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

#include <algorithm>
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

    // Small square, counterclockwise as a standalone simple polygon. Each component
    // (outer boundary or hole) of an MCSCPolygonalDomain is independently normalized to
    // this orientation (mirroring polygon.m's per-polygon calc_angles); the MATLAB
    // reference's outer/inner sign asymmetry is applied later, in intpolys.m/extpolys.m,
    // and is out of scope for this data-only type.
    std::vector<Complex> smallSquareCcw(Complex center, double halfSide)
    {
        return {
            center + Complex(-halfSide, -halfSide),
            center + Complex(halfSide, -halfSide),
            center + Complex(halfSide, halfSide),
            center + Complex(-halfSide, halfSide),
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

    // calc_angles should detect the wrong (clockwise) orientation and flip vertices +
    // angles so the stored data is always counterclockwise, per-component.
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
    // Note: for any closed polygonal path with non-degenerate edges, sum(1-alpha) is
    // always an integer (a topological winding-number quantization) regardless of
    // vertex order or self-intersection -- e.g. a self-intersecting "bowtie"
    // quadrilateral still sums to 0. The only way to break the invariant is a
    // degenerate (zero-length) edge, which makes the turning-angle computation at
    // that vertex ill-defined (angle of a zero vector is taken as 0).
    //
    // A repeated vertex produces exactly such a zero-length edge.
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
    auto hole = smallSquareCcw(Complex(2.0, 2.0), 0.5);

    MCSCPolygonalDomain domain({outer, hole});

    EXPECT_EQ(2, domain.getConnectivity());
    EXPECT_EQ(4, domain.vertexCount(0));
    EXPECT_EQ(4, domain.vertexCount(1));

    ASSERT_EQ(2u, domain.getVertices().size());
    ASSERT_EQ(2u, domain.getAlpha().size());
    EXPECT_EQ(4u, domain.getVertices()[1].size());
    EXPECT_EQ(4u, domain.getAlpha()[1].size());
}

TEST(MCSCPolygonalDomainTest, ThreeComponentConstruction)
{
    auto outer = unitSquareCcw();
    for (auto& v : outer)
    {
        v *= 10.0;
    }
    auto holeA = smallSquareCcw(Complex(2.0, 2.0), 0.5);
    auto holeB = smallSquareCcw(Complex(7.0, 7.0), 0.5);

    MCSCPolygonalDomain domain({outer, holeA, holeB});

    EXPECT_EQ(3, domain.getConnectivity());
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_EQ(4, domain.vertexCount(i));
        EXPECT_EQ(4u, domain.getAlpha(i).size());
    }

    EXPECT_TRUE(domain.contains(Complex(5.0, 5.0)));   // Inside outer, outside both holes.
    EXPECT_FALSE(domain.contains(Complex(2.0, 2.0)));  // Inside holeA.
    EXPECT_FALSE(domain.contains(Complex(7.0, 7.0)));  // Inside holeB.
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

TEST(MCSCPolygonalDomainTest, ContainsUnboundedMultipleComponents)
{
    // Unbounded target domain: the region is the exterior of all disjoint polygons, with
    // no enclosing "outer" component -- a point strictly inside any one polygon must be
    // reported as outside the domain, not just points inside component 0.
    auto squareA = unitSquareCcw();
    auto squareB = smallSquareCcw(Complex(5.0, 5.0), 0.5);

    MCSCPolygonalDomain domain({squareA, squareB}, /*isUnboundedDomain=*/true);

    EXPECT_FALSE(domain.contains(Complex(0.5, 0.5)));  // Inside component 0.
    EXPECT_FALSE(domain.contains(Complex(5.0, 5.0)));  // Inside component 1.
    EXPECT_TRUE(domain.contains(Complex(2.5, 2.5)));   // Outside both.
}

TEST(MCSCPolygonalDomainTest, ContainsRespectsHole)
{
    auto outer = unitSquareCcw();
    for (auto& v : outer)
    {
        v *= 4.0;
    }
    auto hole = smallSquareCcw(Complex(2.0, 2.0), 0.5);

    MCSCPolygonalDomain domain({outer, hole});

    EXPECT_TRUE(domain.contains(Complex(0.5, 0.5)));   // Inside outer, outside hole.
    EXPECT_FALSE(domain.contains(Complex(2.0, 2.0)));  // Inside hole.
    EXPECT_FALSE(domain.contains(Complex(10.0, 10.0))); // Outside outer boundary entirely.
}

TEST(MCSCPolygonalDomainTest, ContainsOnBoundaryPoints)
{
    // Domain::calculateWindingNumber reports points on (or very near) a boundary edge
    // via a sentinel value; MCSCPolygonalDomain::contains must handle this for both the
    // outer boundary and a hole boundary, for bounded and unbounded domains.
    auto outer = unitSquareCcw();
    for (auto& v : outer)
    {
        v *= 4.0;
    }
    auto hole = smallSquareCcw(Complex(2.0, 2.0), 0.5);

    MCSCPolygonalDomain bounded({outer, hole});
    EXPECT_TRUE(bounded.contains(Complex(0.0, 2.0)));   // On the outer boundary.
    EXPECT_FALSE(bounded.contains(Complex(1.5, 2.0)));  // On the hole boundary.

    MCSCPolygonalDomain unbounded({outer}, /*isUnboundedDomain=*/true);
    EXPECT_FALSE(unbounded.contains(Complex(0.0, 2.0))); // On the (only) boundary.
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

TEST(MCSCPolygonalDomainTest, TransformBoundaryFailureLeavesDomainUnmodified)
{
    // transformBoundary must be atomic: a failing transform should leave the domain
    // exactly as it was before the call, not partially applied.
    auto original = unitSquareCcw();
    MCSCPolygonalDomain domain({original});
    auto originalAlpha = domain.getAlpha(0);

    EXPECT_THROW(
        domain.transformBoundary(
            [](const Complex& z) { return (z == Complex(0.0, 0.0)) ? Complex(1.0, 0.0) : z; }),
        std::invalid_argument);

    const auto& vertices = domain.getVertices(0);
    ASSERT_EQ(original.size(), vertices.size());
    for (size_t i = 0; i < original.size(); ++i)
    {
        EXPECT_EQ(original[i], vertices[i]);
    }

    const auto& alpha = domain.getAlpha(0);
    ASSERT_EQ(originalAlpha.size(), alpha.size());
    for (size_t i = 0; i < originalAlpha.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(originalAlpha[i], alpha[i]);
    }
}

TEST(MCSCPolygonalDomainTest, TransformBoundaryReorientingReflection)
{
    // An orientation-reversing transform (e.g. conjugation) does not throw -- it is
    // re-normalized back to counterclockwise, the same as construction-time handling.
    MCSCPolygonalDomain domain({unitSquareCcw()});

    domain.transformBoundary([](const Complex& z) { return std::conj(z); });

    // As with construction-time re-orientation, the flip reorders vertices/angles but
    // does not recompute angle values -- the conjugated (now CW-as-computed, then
    // flipped) square still carries raw turning angles of 1.5.
    const auto& alpha = domain.getAlpha(0);
    ASSERT_EQ(4u, alpha.size());
    for (double a : alpha)
    {
        EXPECT_NEAR(1.5, a, 1e-12);
    }
}

TEST(MCSCPolygonalDomainTest, OutOfRangeComponentIndexThrows)
{
    MCSCPolygonalDomain domain({unitSquareCcw()});

    EXPECT_THROW(domain.getVertices(1), std::invalid_argument);
    EXPECT_THROW(domain.getAlpha(-1), std::invalid_argument);
    EXPECT_THROW(domain.vertexCount(5), std::invalid_argument);
}
