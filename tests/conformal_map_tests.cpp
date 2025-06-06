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
#include "../src/Domain.h"

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

    bool wasComputeCalled() const
    {
        return compute_called;
    }

    double getLastAccuracy() const
    {
        return last_accuracy;
    }
};

// ConformalMap is now concrete, no need for MockMap

TEST(ConformalMapTest, Construction)
{
    // Create domains for the test
    auto source_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);
    auto target_domain = std::make_shared<EllipticalDomain>(
        2.0, 1.0, 0.0, ComplexDouble(0.0, 0.0), false);

    // Test constructor with domains
    ConformalMap map(source_domain, target_domain, nullptr, false);
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
    auto method = std::make_shared<MockMethod>();
    ConformalMap map(source_domain, target_domain, method);

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

    // Create identity map with mock method
    auto method = std::make_shared<MockMethod>();
    ConformalMap map(source_domain, target_domain, method);

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
    ConformalMap map(source_domain, target_domain);

    // Test compute
    method.compute(map, 1e-10);
    EXPECT_EQ(5, method.getIterationCount());
    EXPECT_EQ(1e-10 * 0.5, method.getAccuracy());
}

TEST(ConformalMapMethodTest, ValidateDomainCompatibility)
{
    class TestMethod : public ConformalMapMethod
    {
    public:
        TestMethod() : ConformalMapMethod() {}

        void compute([[maybe_unused]] ConformalMap& map, [[maybe_unused]] double target_accuracy) override
        {
            // Just a stub
        }

        ComplexDouble map(const ComplexDouble& z) const override
        {
            return z;
        }

        ComplexDouble inverseMap(const ComplexDouble& w) const override
        {
            return w;
        }

        void testDomainValidation(std::shared_ptr<Domain> domain, int expected_connectivity)
        {
            validateDomainCompatibility(domain, expected_connectivity);
        }
    };

    TestMethod method;

    // Create domains
    auto simply_connected_domain = std::make_shared<CircularDomain>(
        ComplexDouble(0.0, 0.0), 1.0, false);

    // Test validation with correct connectivity (1 = simply connected)
    EXPECT_NO_THROW(method.testDomainValidation(simply_connected_domain, 1));

    // Test validation with incorrect connectivity - should throw
    EXPECT_THROW(method.testDomainValidation(simply_connected_domain, 0), std::invalid_argument);

    // Test null domain - should throw
    EXPECT_THROW(method.testDomainValidation(nullptr, 0), std::invalid_argument);
}
