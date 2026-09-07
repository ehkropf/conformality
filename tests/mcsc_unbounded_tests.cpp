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

#include "../src/methods/MCSCUnbounded.h"

#include "../src/core/ConformalMap.h"

#include <gtest/gtest.h>
#include <cmath>

namespace
{

// example_driver_ext.m's 3-polygon unbounded example, matching
// mcsc_unbounded_objective_function_tests.cpp's fixture.
MCSCPolygonalDomain makeExampleDriverTargetPolygon()
{
    std::vector<Complex> triangle0 = {Complex(0.0, 0.0), Complex(-1.5, -2.0), Complex(1.5, -2.0)};
    std::vector<Complex> triangle1 = {
        Complex(-1.0, 1.0), Complex(-1.5, 1.0 + std::sqrt(3.0) / 2.0), Complex(-2.0, 1.0)};
    std::vector<Complex> triangle2 = {
        Complex(1.0, 1.0), Complex(2.0, 1.0), Complex(1.5, 1.0 + std::sqrt(3.0) / 2.0)};
    return MCSCPolygonalDomain({triangle0, triangle1, triangle2}, /*isUnboundedDomain=*/true);
}

MCSCCircleDomain makeExampleDriverInitialGuess()
{
    const double twoPiThird = 2.0 * M_PI / 3.0;
    return MCSCCircleDomain(std::vector<MCSCCircleDomain::CircleData>{
        {Complex(0.0, 0.0), 1.0, {0.0, twoPiThird, 2.0 * twoPiThird}},
        {Complex(2.4, 1.3), 0.37, {3.0 * M_PI / 2.0, M_PI / 4.0 + 2.0 * M_PI, 3.0 * M_PI / 4.0 + 2.0 * M_PI}},
        {Complex(2.4, -1.3), 0.37, {M_PI / 2.0, 5.0 * M_PI / 4.0, 7.0 * M_PI / 4.0}},
    });
}

std::shared_ptr<ConformalMap> makeExampleDriverMap()
{
    auto polygon = std::make_shared<MCSCPolygonalDomain>(makeExampleDriverTargetPolygon());
    auto circle = std::make_shared<MCSCCircleDomain>(makeExampleDriverInitialGuess());
    auto method = std::make_shared<MCSCUnbounded>();
    return std::make_shared<ConformalMap>(circle, polygon, method);
}

}  // namespace

TEST(MCSCUnboundedTest, MapThrowsBeforeCompute)
{
    MCSCUnbounded method;
    EXPECT_THROW(method.map(Complex(0.5, 0.0)), std::runtime_error);
}

TEST(MCSCUnboundedTest, InverseMapThrows)
{
    MCSCUnbounded method;
    EXPECT_THROW(method.inverseMap(Complex(0.0, 0.0)), std::runtime_error);
}

TEST(MCSCUnboundedTest, ValidateSourceDomainRejectsWrongType)
{
    MCSCUnbounded method;
    auto polygon = std::make_shared<MCSCPolygonalDomain>(makeExampleDriverTargetPolygon());
    EXPECT_THROW(method.validateDomain(polygon, 3, ConformalMapMethod::DomainRole::Source), std::invalid_argument);
}

TEST(MCSCUnboundedTest, ValidateTargetDomainRejectsWrongType)
{
    MCSCUnbounded method;
    auto circle = std::make_shared<MCSCCircleDomain>(makeExampleDriverInitialGuess());
    EXPECT_THROW(method.validateDomain(circle, 3, ConformalMapMethod::DomainRole::Target), std::invalid_argument);
}

TEST(MCSCUnboundedTest, ComputeSolvesAndVertexImagesMatchTargetPolygon)
{
    auto map = makeExampleDriverMap();
    map->compute();

    auto method = std::dynamic_pointer_cast<MCSCUnbounded>(map->getMethod());
    ASSERT_NE(method, nullptr);

    const MCSCCircleDomain& solved = method->solvedCircleDomain();
    auto polygon = std::dynamic_pointer_cast<MCSCPolygonalDomain>(map->getTargetDomain());
    ASSERT_NE(polygon, nullptr);

    // Evaluating the map at each prevertex should recover the corresponding target vertex --
    // the arc segment is skipped entirely (mindt < 100*eps case) and the result is just the
    // known vertex image, so this is close to a tautology for k=0 but genuinely exercises the
    // full arc-path machinery for k >= 1 on each circle.
    for (int j = 0; j < solved.circleCount(); ++j)
    {
        const auto prevertices = solved.getPrevertices(j);
        const auto& targetVertices = polygon->getVertices(j);
        for (std::size_t k = 0; k < prevertices.size(); ++k)
        {
            const Complex w = map->map(prevertices[k]);
            EXPECT_NEAR(w.real(), targetVertices[k].real(), 1e-6);
            EXPECT_NEAR(w.imag(), targetVertices[k].imag(), 1e-6);
        }
    }
}

TEST(MCSCUnboundedTest, MapAtBoundaryPointSkipsLineSegment)
{
    auto map = makeExampleDriverMap();
    map->compute();

    auto method = std::dynamic_pointer_cast<MCSCUnbounded>(map->getMethod());
    ASSERT_NE(method, nullptr);
    const MCSCCircleDomain& solved = method->solvedCircleDomain();

    // A point on circle 0's boundary, midway between two prevertices in angle -- exercises the
    // arc-only path (line segment skipped).
    const auto& t = solved.getPrevertexAngles(0);
    const double midAngle = 0.5 * (t[0] + t[1]);
    const Complex z = solved.getCenter(0) + solved.getRadius(0) * std::exp(Complex(0.0, 1.0) * midAngle);

    EXPECT_NO_THROW({
        Complex w = map->map(z);
        EXPECT_TRUE(std::isfinite(w.real()));
        EXPECT_TRUE(std::isfinite(w.imag()));
    });
}

TEST(MCSCUnboundedTest, MapIsDeterministic)
{
    auto map = makeExampleDriverMap();
    map->compute();

    const Complex z(3.0, 0.5);
    const Complex w1 = map->map(z);
    const Complex w2 = map->map(z);

    EXPECT_DOUBLE_EQ(w1.real(), w2.real());
    EXPECT_DOUBLE_EQ(w1.imag(), w2.imag());
}

TEST(MCSCUnboundedTest, MapAtGeneralExteriorPointIsFinite)
{
    auto map = makeExampleDriverMap();
    map->compute();

    const Complex z(4.0, 3.0);
    const Complex w = map->map(z);

    EXPECT_TRUE(std::isfinite(w.real()));
    EXPECT_TRUE(std::isfinite(w.imag()));
}
