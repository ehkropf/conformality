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
#include <memory>
#include "../src/ConformalMap.h"
#include "../src/ConformalMapMethod.h"
#include "../src/Domain.hpp"

// Mock implementation of ConformalMapMethod for testing
class MockMethod : public ConformalMapMethod
{
private:
    bool compute_called = false;
    double last_accuracy = 0.0;

public:
    MockMethod() : ConformalMapMethod() {}

    void compute([[maybe_unused]] ConformalMap& map, double target_accuracy) override
    {
        compute_called = true;
        last_accuracy = target_accuracy;
        iteration_count = 5; // Simulate some iterations
        accuracy = target_accuracy * 0.5; // Simulate achieved accuracy
    }

    bool wasComputeCalled() const
    {
        return compute_called;
    }

    double getLastAccuracy() const
    {
        return last_accuracy;
    }
};

// Mock implementation of ConformalMap for testing
class MockMap : public ConformalMap
{
public:
    MockMap(std::shared_ptr<Domain> source, std::shared_ptr<Domain> target, bool external = false)
        : ConformalMap(source, target, external)
    {
    }

    ComplexDouble map(const ComplexDouble& z) const override
    {
        // Simple identity mapping for testing
        return z;
    }

    ComplexDouble inverseMap(const ComplexDouble& w) const override
    {
        // Simple identity mapping for testing
        return w;
    }
};

TEST(ConformalMapTest, Construction)
{
    // Create domains for the test
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<EllipticalDomain>(
        2.0, 1.0, 0.0, ComplexDouble(0.0, 0.0), false);

    // Test constructor with domains
    MockMap map(source_domain, target_domain, false);
    EXPECT_EQ(source_domain, map.getSourceDomain());
    EXPECT_EQ(target_domain, map.getTargetDomain());
    EXPECT_FALSE(map.isExternalMap());

    // Test setting external flag
    map.setExternal(true);
    EXPECT_TRUE(map.isExternalMap());
}

TEST(ConformalMapTest, MethodInteraction)
{
    // Create domains for the test
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<EllipticalDomain>(
        2.0, 1.0, 0.0, ComplexDouble(0.0, 0.0), false);

    // Create map and method
    MockMap map(source_domain, target_domain);
    auto method = std::make_shared<MockMethod>();

    // Test setting method
    map.setMethod(method);
    EXPECT_EQ(method, map.getMethod());

    // Test compute delegation
    map.compute(1e-8);
    EXPECT_TRUE(std::dynamic_pointer_cast<MockMethod>(map.getMethod())->wasComputeCalled());
    EXPECT_EQ(1e-8, std::dynamic_pointer_cast<MockMethod>(map.getMethod())->getLastAccuracy());
}

TEST(ConformalMapTest, Identity)
{
    // Create domains for the test
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);

    // Create identity map
    MockMap map(source_domain, target_domain);

    // Test identity mapping
    ComplexDouble z(0.5, 0.5);
    ComplexDouble w = map.map(z);
    EXPECT_EQ(z, w);

    ComplexDouble z_back = map.inverseMap(w);
    EXPECT_EQ(z, z_back);
}

TEST(ConformalMapMethodTest, BasicProperties)
{
    MockMethod method;

    // Test default values
    EXPECT_EQ(0.0, method.getAccuracy());
    EXPECT_EQ(0, method.getIterationCount());

    // Test setting max iterations
    method.setMaxIterations(100);

    // Create domains and map for testing compute
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<EllipticalDomain>(
        2.0, 1.0, 0.0, ComplexDouble(0.0, 0.0), false);
    MockMap map(source_domain, target_domain);

    // Test compute
    method.compute(map, 1e-10);
    EXPECT_EQ(5, method.getIterationCount());
    EXPECT_EQ(1e-10 * 0.5, method.getAccuracy());
}

TEST(ConformalMapMethodTest, ValidateMapType)
{
    class TestMethod : public ConformalMapMethod
    {
    public:
        TestMethod() : ConformalMapMethod() {}

        void compute([[maybe_unused]] ConformalMap& map, [[maybe_unused]] double target_accuracy) override
        {
            // Just a stub
        }

        void testValidation(ConformalMap& map, const std::string& expected_type)
        {
            validateMapType(map, expected_type);
        }
    };

    TestMethod method;

    // Create domains and map
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<EllipticalDomain>(
        2.0, 1.0, 0.0, ComplexDouble(0.0, 0.0), false);
    MockMap map(source_domain, target_domain);

    // Test validation with correct type
    EXPECT_NO_THROW(method.testValidation(map, "MockMap"));

    // Test validation with incorrect type - should throw
    EXPECT_THROW(method.testValidation(map, "WrongMapType"), std::invalid_argument);
}
