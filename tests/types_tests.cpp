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
#include "../src/Domain.h"

class TypesTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Set up common test data
        bounded_source = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 1.0, false);
        bounded_target = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 2.0, false);
        unbounded_source = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 1.0, true);
        unbounded_target = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 2.0, true);
    }

    std::shared_ptr<CircularDomain> bounded_source;
    std::shared_ptr<CircularDomain> bounded_target;
    std::shared_ptr<CircularDomain> unbounded_source;
    std::shared_ptr<CircularDomain> unbounded_target;
};

TEST_F(TypesTest, MappingTypeEnum)
{
    // Test that enum values are as expected
    EXPECT_EQ(static_cast<int>(MappingType::INTERIOR_TO_INTERIOR), 0);
    EXPECT_EQ(static_cast<int>(MappingType::EXTERIOR_TO_INTERIOR), 1);
    EXPECT_EQ(static_cast<int>(MappingType::INTERIOR_TO_EXTERIOR), 2);
    EXPECT_EQ(static_cast<int>(MappingType::EXTERIOR_TO_EXTERIOR), 3);
}

TEST_F(TypesTest, DetermineMappingTypeBoundedToBounded)
{
    MappingType result = determineMappingType(*bounded_source, *bounded_target);
    EXPECT_EQ(result, MappingType::INTERIOR_TO_INTERIOR);
}

TEST_F(TypesTest, DetermineMappingTypeUnboundedToBounded)
{
    MappingType result = determineMappingType(*unbounded_source, *bounded_target);
    EXPECT_EQ(result, MappingType::EXTERIOR_TO_INTERIOR);
}

TEST_F(TypesTest, DetermineMappingTypeBoundedToUnbounded)
{
    MappingType result = determineMappingType(*bounded_source, *unbounded_target);
    EXPECT_EQ(result, MappingType::INTERIOR_TO_EXTERIOR);
}

TEST_F(TypesTest, DetermineMappingTypeUnboundedToUnbounded)
{
    MappingType result = determineMappingType(*unbounded_source, *unbounded_target);
    EXPECT_EQ(result, MappingType::EXTERIOR_TO_EXTERIOR);
}

TEST_F(TypesTest, DomainUnboundedProperty)
{
    // Test that domains report their boundedness correctly
    EXPECT_FALSE(bounded_source->isUnbounded());
    EXPECT_FALSE(bounded_target->isUnbounded());
    EXPECT_TRUE(unbounded_source->isUnbounded());
    EXPECT_TRUE(unbounded_target->isUnbounded());
}

TEST_F(TypesTest, MappingTypeWithDifferentDomainTypes)
{
    // Test with elliptical domains
    auto bounded_ellipse = std::make_shared<EllipticalDomain>(2.0, 1.0, 0.0, Complex(0.0, 0.0), false);
    auto unbounded_ellipse = std::make_shared<EllipticalDomain>(2.0, 1.0, 0.0, Complex(0.0, 0.0), true);
    
    // Mix different domain types
    EXPECT_EQ(determineMappingType(*bounded_ellipse, *bounded_source), MappingType::INTERIOR_TO_INTERIOR);
    EXPECT_EQ(determineMappingType(*unbounded_ellipse, *bounded_source), MappingType::EXTERIOR_TO_INTERIOR);
    EXPECT_EQ(determineMappingType(*bounded_ellipse, *unbounded_source), MappingType::INTERIOR_TO_EXTERIOR);
    EXPECT_EQ(determineMappingType(*unbounded_ellipse, *unbounded_source), MappingType::EXTERIOR_TO_EXTERIOR);
}

TEST_F(TypesTest, ComplexTypeAlias)
{
    // Test that Complex type alias works correctly
    Complex z1(1.0, 2.0);
    Complex z2 = std::complex<double>(3.0, 4.0);
    
    EXPECT_DOUBLE_EQ(std::real(z1), 1.0);
    EXPECT_DOUBLE_EQ(std::imag(z1), 2.0);
    EXPECT_DOUBLE_EQ(std::real(z2), 3.0);
    EXPECT_DOUBLE_EQ(std::imag(z2), 4.0);
    
    Complex sum = z1 + z2;
    EXPECT_DOUBLE_EQ(std::real(sum), 4.0);
    EXPECT_DOUBLE_EQ(std::imag(sum), 6.0);
}

TEST_F(TypesTest, BoundaryToleranceConstant)
{
    // Test that the boundary tolerance constant is defined and reasonable
    EXPECT_GT(BOUNDARY_TOLERANCE, 0.0);
    EXPECT_LT(BOUNDARY_TOLERANCE, 1e-6); // Should be very small
    EXPECT_EQ(BOUNDARY_TOLERANCE, 1e-12); // Should be exactly what we expect
}