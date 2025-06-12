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
#include "../src/Domain.h"
#include "../src/Boundary.h"
#include "../src/BoundaryComponent.h"
#include <memory>
#include <cmath>

class SimplyConnectedDomainTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a unit circle boundary for testing
        auto circleFunc = [](double t) -> Complex
        {
            return Complex(std::cos(t), std::sin(t));
        };
        auto circleDerivFunc = [](double t) -> Complex
        {
            return Complex(-std::sin(t), std::cos(t));
        };

        circleComponent = std::make_shared<AnalyticBoundaryComponent>(circleFunc, circleDerivFunc);
        circleBoundary = std::make_shared<Boundary>(circleComponent);

        // Create an ellipse boundary for testing
        double a = 2.0, b = 1.0;
        auto ellipseFunc = [a, b](double t) -> Complex
        {
            return Complex(a * std::cos(t), b * std::sin(t));
        };
        auto ellipseDerivFunc = [a, b](double t) -> Complex
        {
            return Complex(-a * std::sin(t), b * std::cos(t));
        };

        ellipseComponent = std::make_shared<AnalyticBoundaryComponent>(ellipseFunc, ellipseDerivFunc);
        ellipseBoundary = std::make_shared<Boundary>(ellipseComponent);

        // Create a square boundary using discrete points
        std::vector<Complex> squarePoints = {
            Complex(1.0, 1.0),
            Complex(-1.0, 1.0),
            Complex(-1.0, -1.0),
            Complex(1.0, -1.0)
        };
        squareComponent = std::make_shared<DiscreteBoundaryComponent>(squarePoints);
        squareBoundary = std::make_shared<Boundary>(squareComponent);
    }

    std::shared_ptr<AnalyticBoundaryComponent> circleComponent{nullptr};
    std::shared_ptr<Boundary> circleBoundary{nullptr};
    std::shared_ptr<AnalyticBoundaryComponent> ellipseComponent{nullptr};
    std::shared_ptr<Boundary> ellipseBoundary{nullptr};
    std::shared_ptr<DiscreteBoundaryComponent> squareComponent{nullptr};
    std::shared_ptr<Boundary> squareBoundary{nullptr};
};

TEST_F(SimplyConnectedDomainTest, Construction)
{
    // Test internal domain construction
    SimplyConnectedDomain internalDomain(circleBoundary, false);

    EXPECT_FALSE(internalDomain.isUnbounded());
    EXPECT_EQ(1, internalDomain.getConnectivity());
    EXPECT_EQ(circleBoundary.get(), &internalDomain.getBoundary());

    // Test external domain construction
    SimplyConnectedDomain externalDomain(circleBoundary, true);

    EXPECT_TRUE(externalDomain.isUnbounded());
    EXPECT_EQ(1, externalDomain.getConnectivity());
    EXPECT_EQ(circleBoundary.get(), &externalDomain.getBoundary());
}

TEST_F(SimplyConnectedDomainTest, BoundaryAccess)
{
    SimplyConnectedDomain domain(ellipseBoundary, false);

    // Test const boundary access
    const SimplyConnectedDomain& constDomain = domain;
    const Boundary& constBoundary = constDomain.getBoundary();
    EXPECT_EQ(ellipseBoundary.get(), &constBoundary);

    // Test non-const boundary access
    Boundary& boundary = domain.getBoundary();
    EXPECT_EQ(ellipseBoundary.get(), &boundary);

    // Verify we can sample from the boundary
    auto samples = boundary.sample(8);
    EXPECT_EQ(1, samples.size()); // One component
    EXPECT_EQ(8, samples[0].size()); // 8 points in that component
}

TEST_F(SimplyConnectedDomainTest, InternalCircularDomainContainment)
{
    SimplyConnectedDomain domain(circleBoundary, false);

    // Test points clearly inside the unit circle
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.0)));     // Origin
    EXPECT_TRUE(domain.contains(Complex(0.5, 0.0)));     // On x-axis
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.5)));     // On y-axis
    EXPECT_TRUE(domain.contains(Complex(0.3, 0.4)));     // General interior point
    EXPECT_TRUE(domain.contains(Complex(-0.2, -0.3)));   // Negative quadrant

    // Test points clearly outside the unit circle
    EXPECT_FALSE(domain.contains(Complex(2.0, 0.0)));    // Far on x-axis
    EXPECT_FALSE(domain.contains(Complex(0.0, 2.0)));    // Far on y-axis
    EXPECT_FALSE(domain.contains(Complex(1.5, 1.5)));    // Far in positive quadrant
    EXPECT_FALSE(domain.contains(Complex(-2.0, -2.0)));  // Far in negative quadrant

    // Test points near the boundary (should be inside for internal domain)
    EXPECT_TRUE(domain.contains(Complex(0.95, 0.0)));    // Close to boundary on x-axis
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.95)));    // Close to boundary on y-axis
}

TEST_F(SimplyConnectedDomainTest, ExternalCircularDomainContainment)
{
    SimplyConnectedDomain domain(circleBoundary, true);

    // Test points clearly inside the unit circle (should be outside external domain)
    EXPECT_FALSE(domain.contains(Complex(0.0, 0.0)));    // Origin
    EXPECT_FALSE(domain.contains(Complex(0.5, 0.0)));    // On x-axis
    EXPECT_FALSE(domain.contains(Complex(0.0, 0.5)));    // On y-axis
    EXPECT_FALSE(domain.contains(Complex(0.3, 0.4)));    // General interior point

    // Test points clearly outside the unit circle (should be inside external domain)
    EXPECT_TRUE(domain.contains(Complex(2.0, 0.0)));     // Far on x-axis
    EXPECT_TRUE(domain.contains(Complex(0.0, 2.0)));     // Far on y-axis
    EXPECT_TRUE(domain.contains(Complex(1.5, 1.5)));     // Far in positive quadrant
    EXPECT_TRUE(domain.contains(Complex(-2.0, -2.0)));   // Far in negative quadrant

    // Test points near the boundary from outside
    EXPECT_TRUE(domain.contains(Complex(1.05, 0.0)));    // Just outside boundary on x-axis
    EXPECT_TRUE(domain.contains(Complex(0.0, 1.05)));    // Just outside boundary on y-axis
}

TEST_F(SimplyConnectedDomainTest, EllipticalDomainContainment)
{
    SimplyConnectedDomain domain(ellipseBoundary, false);

    // Test points clearly inside the ellipse (a=2, b=1)
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.0)));     // Origin
    EXPECT_TRUE(domain.contains(Complex(1.0, 0.0)));     // Inside along major axis
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.5)));     // Inside along minor axis
    EXPECT_TRUE(domain.contains(Complex(0.5, 0.3)));     // General interior point

    // Test points clearly outside the ellipse
    EXPECT_FALSE(domain.contains(Complex(3.0, 0.0)));    // Outside along major axis
    EXPECT_FALSE(domain.contains(Complex(0.0, 2.0)));    // Outside along minor axis
    EXPECT_FALSE(domain.contains(Complex(2.0, 1.0)));    // Outside in general position

    // Test points near the boundary
    EXPECT_TRUE(domain.contains(Complex(1.8, 0.0)));     // Near major axis endpoint
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.9)));     // Near minor axis endpoint
}

TEST_F(SimplyConnectedDomainTest, DiscreteBoundaryContainment)
{
    SimplyConnectedDomain domain(squareBoundary, false);

    // Test points clearly inside the square
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.0)));     // Origin
    EXPECT_TRUE(domain.contains(Complex(0.5, 0.5)));     // Interior point
    EXPECT_TRUE(domain.contains(Complex(-0.5, 0.5)));    // Another interior point

    // Test points clearly outside the square
    EXPECT_FALSE(domain.contains(Complex(2.0, 0.0)));    // Outside right
    EXPECT_FALSE(domain.contains(Complex(0.0, 2.0)));    // Outside top
    EXPECT_FALSE(domain.contains(Complex(-2.0, -2.0)));  // Outside bottom-left

    // Test points near the boundary
    EXPECT_TRUE(domain.contains(Complex(0.9, 0.9)));     // Close to corner
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.9)));     // Close to edge
}

TEST_F(SimplyConnectedDomainTest, BoundaryTransformation)
{
    SimplyConnectedDomain domain(circleBoundary, false);

    // Define a simple transformation (scaling by factor of 2)
    auto scaleTransform = [](const Complex& z) -> Complex
    {
        return z * Complex(2.0, 0.0);
    };

    // Get some points on the original boundary
    auto originalSamplesVec = domain.getBoundary().sample(8);
    std::vector<Complex> originalSamples = originalSamplesVec[0]; // Single component

    // Verify original points are on unit circle
    for (const Complex& point : originalSamples)
    {
        EXPECT_NEAR(1.0, std::abs(point), 1e-10);
    }

    // Apply transformation
    domain.transformBoundary(scaleTransform);

    // Get points on the transformed boundary
    auto transformedSamplesVec = domain.getBoundary().sample(8);
    std::vector<Complex> transformedSamples = transformedSamplesVec[0]; // Single component

    // Verify transformed points are on circle of radius 2
    for (const Complex& point : transformedSamples)
    {
        EXPECT_NEAR(2.0, std::abs(point), 1e-1); // Relaxed tolerance due to discrete approximation
    }

    // Test containment after transformation
    EXPECT_TRUE(domain.contains(Complex(1.5, 0.0)));     // Inside scaled circle
    EXPECT_TRUE(domain.contains(Complex(0.5, 0.0)));     // Also inside
}

TEST_F(SimplyConnectedDomainTest, ComplexTransformation)
{
    SimplyConnectedDomain domain(circleBoundary, false);

    // Define a more complex transformation (z -> z^2)
    auto squareTransform = [](const Complex& z) -> Complex
    {
        return z * z;
    };

    // Apply transformation
    domain.transformBoundary(squareTransform);

    // After z -> z^2 transformation, the unit circle should become a more complex shape
    // Test some expected behavior
    auto transformedSamplesVec = domain.getBoundary().sample(16);
    EXPECT_EQ(1, transformedSamplesVec.size()); // One component
    auto transformedSamples = transformedSamplesVec[0];
    EXPECT_EQ(16, transformedSamples.size());

    // The transformation should preserve some symmetries
    // Points should still exist, though the shape will be different
    bool hasPositiveRealPoints = false;
    bool hasNegativeRealPoints = false;
    bool hasPositiveImagPoints = false;
    bool hasNegativeImagPoints = false;

    for (const auto& point : transformedSamples)
    {
        if (point.real() > 0.1) hasPositiveRealPoints = true;
        if (point.real() < -0.1) hasNegativeRealPoints = true;
        if (point.imag() > 0.1) hasPositiveImagPoints = true;
        if (point.imag() < -0.1) hasNegativeImagPoints = true;
    }

    // After z^2, we should have points in multiple quadrants
    EXPECT_TRUE(hasPositiveRealPoints);
    EXPECT_TRUE(hasNegativeRealPoints);
    EXPECT_TRUE(hasPositiveImagPoints);
    EXPECT_TRUE(hasNegativeImagPoints);
}

TEST_F(SimplyConnectedDomainTest, EdgeCasesAndRobustness)
{
    SimplyConnectedDomain domain(circleBoundary, false);

    // Test with points exactly on the boundary (should be considered inside for numerical stability)
    EXPECT_TRUE(domain.contains(Complex(1.0, 0.0)));     // Exactly on boundary
    EXPECT_TRUE(domain.contains(Complex(0.0, 1.0)));     // Exactly on boundary
    EXPECT_TRUE(domain.contains(Complex(-1.0, 0.0)));    // Exactly on boundary

    // Test boundary access after construction
    const Boundary& boundary = domain.getBoundary();
    EXPECT_EQ(1, boundary.getNumComponents());

    const BoundaryComponent& component = boundary.getComponent(0);
    Complex testPoint = component.evaluate(0.0);
    EXPECT_NEAR(1.0, testPoint.real(), 1e-10);
    EXPECT_NEAR(0.0, testPoint.imag(), 1e-10);
}

TEST_F(SimplyConnectedDomainTest, IdentityTransformation)
{
    SimplyConnectedDomain domain(ellipseBoundary, false);

    // Store original state
    auto originalSamplesVec = domain.getBoundary().sample(10);
    auto originalSamples = originalSamplesVec[0]; // Single component

    // Apply identity transformation
    auto identityTransform = [](const Complex& z) -> Complex
    {
        return z;
    };

    domain.transformBoundary(identityTransform);

    // Check that domain behavior is preserved
    auto newSamplesVec = domain.getBoundary().sample(10);
    auto newSamples = newSamplesVec[0]; // Single component
    EXPECT_EQ(originalSamples.size(), newSamples.size());

    // The discrete representation might not be exactly identical due to resampling,
    // but the containment behavior should be similar
    EXPECT_TRUE(domain.contains(Complex(0.0, 0.0)));     // Still inside
    EXPECT_FALSE(domain.contains(Complex(3.0, 0.0)));    // Still outside
}
