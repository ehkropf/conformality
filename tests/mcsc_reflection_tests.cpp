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

#include "../src/methods/MCSCReflection.h"

#include <gtest/gtest.h>
#include <cmath>

using mcsc::reflectCircle;
using mcsc::reflectCircleSequence;
using mcsc::reflectPoint;
using mcsc::ReflectedCircle;

TEST(MCSCReflectionTest, ReflectPointThroughUnitCircleAtOrigin)
{
    // z = 2 reflected through the unit circle at 0 should land at 0.5 (real axis, standard
    // circle inversion).
    Complex zo = reflectPoint(Complex(0.0, 0.0), 1.0, Complex(2.0, 0.0));
    EXPECT_NEAR(std::real(zo), 0.5, 1e-12);
    EXPECT_NEAR(std::imag(zo), 0.0, 1e-12);
}

TEST(MCSCReflectionTest, ReflectPointIsSelfInverse)
{
    // Reflecting twice through the same circle returns the original point.
    Complex c(0.3, -0.2);
    double r = 0.7;
    Complex z(1.5, 0.4);

    Complex once = reflectPoint(c, r, z);
    Complex twice = reflectPoint(c, r, once);

    EXPECT_NEAR(std::real(twice), std::real(z), 1e-10);
    EXPECT_NEAR(std::imag(twice), std::imag(z), 1e-10);
}

TEST(MCSCReflectionTest, ReflectPointThrowsOnNonPositiveRadius)
{
    EXPECT_THROW(reflectPoint(Complex(0.0, 0.0), 0.0, Complex(1.0, 0.0)), std::invalid_argument);
    EXPECT_THROW(reflectPoint(Complex(0.0, 0.0), -1.0, Complex(1.0, 0.0)), std::invalid_argument);
}

TEST(MCSCReflectionTest, ReflectPointThrowsAtReflectingCenter)
{
    EXPECT_THROW(reflectPoint(Complex(0.0, 0.0), 1.0, Complex(0.0, 0.0)), std::invalid_argument);
}

TEST(MCSCReflectionTest, ReflectCircleExternalToUnitCircleAtOrigin)
{
    // Circle centered at 3, radius 0.5, reflected through the unit circle at 0.
    // co = 0 + 1*(3-0)/(9-0.25) = 3/8.75; ro = 1*0.5/8.75.
    auto [co, ro] = reflectCircle(Complex(0.0, 0.0), 1.0, Complex(3.0, 0.0), 0.5);
    EXPECT_NEAR(std::real(co), 3.0 / 8.75, 1e-12);
    EXPECT_NEAR(std::imag(co), 0.0, 1e-12);
    EXPECT_NEAR(ro, 0.5 / 8.75, 1e-12);
}

TEST(MCSCReflectionTest, ReflectCircleThrowsOnNonPositiveRadius)
{
    EXPECT_THROW(reflectCircle(Complex(0.0, 0.0), 0.0, Complex(3.0, 0.0), 0.5), std::invalid_argument);
    EXPECT_THROW(reflectCircle(Complex(0.0, 0.0), 1.0, Complex(3.0, 0.0), 0.0), std::invalid_argument);
}

TEST(MCSCReflectionTest, ReflectCircleHandlesCircleContainingReflectingCircle)
{
    // (ci, ri) properly containing (c, r) -- |ci - c| < ri -- makes denom = |ci-c|^2 - ri^2
    // negative, a legitimate configuration the abs(denom) in the radius formula must still
    // handle correctly (as opposed to only ever exercising denom > 0, the more common case
    // above where (ci, ri) is exterior to (c, r)).
    // co = 0 + 1*(1-0)/(1-4) = -1/3; ro = 1*2/3 = 2/3.
    auto [co, ro] = reflectCircle(Complex(0.0, 0.0), 1.0, Complex(1.0, 0.0), 2.0);
    EXPECT_NEAR(std::real(co), -1.0 / 3.0, 1e-12);
    EXPECT_NEAR(std::imag(co), 0.0, 1e-12);
    EXPECT_NEAR(ro, 2.0 / 3.0, 1e-12);
}

TEST(MCSCReflectionTest, SequenceLevelZeroIsOriginalCircles)
{
    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0), Complex(-3.0, 0.0)};
    std::vector<double> radii = {1.0, 0.5, 0.5};
    std::vector<std::vector<Complex>> prevertices = {
        {Complex(1.0, 0.0), Complex(-1.0, 0.0)},
        {Complex(3.5, 0.0), Complex(2.5, 0.0)},
        {Complex(-2.5, 0.0), Complex(-3.5, 0.0)},
    };
    std::vector<Complex> outerCenters = centers;

    auto result = reflectCircleSequence(centers, radii, prevertices, outerCenters, 0);

    ASSERT_EQ(result.size(), 3u);
    for (std::size_t j = 0; j < 3; ++j)
    {
        ASSERT_EQ(result[j].size(), 1u);
        EXPECT_EQ(result[j][0].center, centers[j]);
        EXPECT_DOUBLE_EQ(result[j][0].radius, radii[j]);
        EXPECT_EQ(result[j][0].prevertices, prevertices[j]);
        EXPECT_EQ(result[j][0].outerImage, outerCenters[j]);
        EXPECT_EQ(result[j][0].lastReflectedThrough, static_cast<int>(j));
    }
}

TEST(MCSCReflectionTest, SequenceLevelOneBranchCountForThreeCircles)
{
    // m = 3 circles: level 0 has 1 entry per circle, level 1 adds (m-1) = 2 children each.
    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0), Complex(-3.0, 0.0)};
    std::vector<double> radii = {1.0, 0.5, 0.5};
    std::vector<std::vector<Complex>> prevertices = {
        {Complex(1.0, 0.0), Complex(-1.0, 0.0)},
        {Complex(3.5, 0.0), Complex(2.5, 0.0)},
        {Complex(-2.5, 0.0), Complex(-3.5, 0.0)},
    };
    std::vector<Complex> outerCenters = centers;

    auto result = reflectCircleSequence(centers, radii, prevertices, outerCenters, 1);

    for (std::size_t j = 0; j < 3; ++j)
    {
        EXPECT_EQ(result[j].size(), 1u + 2u);  // level 0 (1) + level 1 (m-1=2)
    }
}

TEST(MCSCReflectionTest, SequenceLevelOneNeverReflectsThroughOwnCircle)
{
    // Level-0 entries have lastReflectedThrough == j by construction (the original,
    // unreflected circle) -- that is expected, not a violation. The actual invariant
    // (reflectzsmi's jlr bookkeeping) is that level >= 1 entries are never reflected through
    // the circle they were just produced from, i.e. level 1 (immediately after the original)
    // is never reflected through circle j itself.
    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0), Complex(-3.0, 0.0)};
    std::vector<double> radii = {1.0, 0.5, 0.5};
    std::vector<std::vector<Complex>> prevertices = {
        {Complex(1.0, 0.0), Complex(-1.0, 0.0)},
        {Complex(3.5, 0.0), Complex(2.5, 0.0)},
        {Complex(-2.5, 0.0), Complex(-3.5, 0.0)},
    };
    std::vector<Complex> outerCenters = centers;

    auto result = reflectCircleSequence(centers, radii, prevertices, outerCenters, 2);

    // Level 1 for circle j is indices [1, m) (m-1 = 2 entries); none should be reflected
    // through circle j itself.
    for (std::size_t j = 0; j < 3; ++j)
    {
        for (std::size_t nu = 1; nu < 3; ++nu)
        {
            EXPECT_NE(result[j][nu].lastReflectedThrough, static_cast<int>(j));
        }
    }
}

TEST(MCSCReflectionTest, SequenceThrowsOnFewerThanTwoCircles)
{
    std::vector<Complex> centers = {Complex(0.0, 0.0)};
    std::vector<double> radii = {1.0};
    std::vector<std::vector<Complex>> prevertices = {{Complex(1.0, 0.0)}};
    std::vector<Complex> outerCenters = centers;

    EXPECT_THROW(reflectCircleSequence(centers, radii, prevertices, outerCenters, 0), std::invalid_argument);
}

TEST(MCSCReflectionTest, SequenceThrowsOnEmptyCircleList)
{
    EXPECT_THROW(
        reflectCircleSequence(
            std::vector<Complex>{}, std::vector<double>{}, std::vector<std::vector<Complex>>{},
            std::vector<Complex>{}, 0),
        std::invalid_argument);
}

TEST(MCSCReflectionTest, SequenceThrowsOnNegativeN)
{
    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0)};
    std::vector<double> radii = {1.0, 0.5};
    std::vector<std::vector<Complex>> prevertices = {{Complex(1.0, 0.0)}, {Complex(3.5, 0.0)}};
    std::vector<Complex> outerCenters = centers;

    EXPECT_THROW(reflectCircleSequence(centers, radii, prevertices, outerCenters, -1), std::invalid_argument);
}

TEST(MCSCReflectionTest, SequenceThrowsOnSizeMismatch)
{
    std::vector<Complex> centers = {Complex(0.0, 0.0), Complex(3.0, 0.0)};
    std::vector<double> radii = {1.0};  // wrong size
    std::vector<std::vector<Complex>> prevertices = {{Complex(1.0, 0.0)}, {Complex(3.5, 0.0)}};
    std::vector<Complex> outerCenters = centers;

    EXPECT_THROW(reflectCircleSequence(centers, radii, prevertices, outerCenters, 0), std::invalid_argument);
}
