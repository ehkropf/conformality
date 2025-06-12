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
#include "../src/Boundary.h"
#include "../src/BoundaryComponent.h"

class BoundaryTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a circle component
        auto circleFunc = [](double t) -> Complex {
            return Complex(std::cos(t), std::sin(t));
        };
        auto circleDerivFunc = [](double t) -> Complex {
            return Complex(-std::sin(t), std::cos(t));
        };
        m_circleComponent = std::make_shared<AnalyticBoundaryComponent>(circleFunc, circleDerivFunc);

        // Create an ellipse component
        double a = 2.0, b = 1.5;
        auto ellipseFunc = [a, b](double t) -> Complex {
            return Complex(a * std::cos(t), b * std::sin(t));
        };
        auto ellipseDerivFunc = [a, b](double t) -> Complex {
            return Complex(-a * std::sin(t), b * std::cos(t));
        };
        m_ellipseComponent = std::make_shared<AnalyticBoundaryComponent>(ellipseFunc, ellipseDerivFunc);
    }

    std::shared_ptr<BoundaryComponent> m_circleComponent;
    std::shared_ptr<BoundaryComponent> m_ellipseComponent;
};

TEST_F(BoundaryTest, DefaultConstructor)
{
    Boundary boundary;
    EXPECT_EQ(0, boundary.getNumComponents());
}

TEST_F(BoundaryTest, SingleComponentConstructor)
{
    Boundary boundary(m_circleComponent);
    EXPECT_EQ(1, boundary.getNumComponents());
}

TEST_F(BoundaryTest, MultipleComponentsConstructor)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);
    EXPECT_EQ(2, boundary.getNumComponents());
}

TEST_F(BoundaryTest, AddComponent)
{
    Boundary boundary;
    boundary.addComponent(m_circleComponent);
    EXPECT_EQ(1, boundary.getNumComponents());

    boundary.addComponent(m_ellipseComponent);
    EXPECT_EQ(2, boundary.getNumComponents());
}

TEST_F(BoundaryTest, EvaluateDefaultComponent)
{
    Boundary boundary(m_circleComponent);

    Complex pt1 = boundary.evaluate(0.0);
    EXPECT_NEAR(1.0, pt1.real(), 1e-10);
    EXPECT_NEAR(0.0, pt1.imag(), 1e-10);

    Complex pt2 = boundary.evaluate(M_PI / 2);
    EXPECT_NEAR(0.0, pt2.real(), 1e-10);
    EXPECT_NEAR(1.0, pt2.imag(), 1e-10);
}

TEST_F(BoundaryTest, EvaluateSpecificComponent)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    // Test first component (circle)
    Complex pt1 = boundary.evaluate(0.0, 0);
    EXPECT_NEAR(1.0, pt1.real(), 1e-10);
    EXPECT_NEAR(0.0, pt1.imag(), 1e-10);

    // Test second component (ellipse with a=2.0, b=1.5)
    Complex pt2 = boundary.evaluate(0.0, 1);
    EXPECT_NEAR(2.0, pt2.real(), 1e-10);
    EXPECT_NEAR(0.0, pt2.imag(), 1e-10);
}

TEST_F(BoundaryTest, EvaluateDerivativeDefaultComponent)
{
    Boundary boundary(m_circleComponent);

    Complex deriv1 = boundary.evaluateDerivative(0.0);
    EXPECT_NEAR(0.0, deriv1.real(), 1e-10);
    EXPECT_NEAR(1.0, deriv1.imag(), 1e-10);
}

TEST_F(BoundaryTest, EvaluateDerivativeSpecificComponent)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    // Test first component (circle)
    Complex deriv1 = boundary.evaluateDerivative(0.0, 0);
    EXPECT_NEAR(0.0, deriv1.real(), 1e-10);
    EXPECT_NEAR(1.0, deriv1.imag(), 1e-10);

    // Test second component (ellipse with a=2.0, b=1.5)
    Complex deriv2 = boundary.evaluateDerivative(0.0, 1);
    EXPECT_NEAR(0.0, deriv2.real(), 1e-10);
    EXPECT_NEAR(1.5, deriv2.imag(), 1e-10);
}

TEST_F(BoundaryTest, Sample)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    size_t numPoints = 4;
    auto samples = boundary.sample(numPoints);

    EXPECT_EQ(2, samples.size()); // Two components
    EXPECT_EQ(numPoints, samples[0].size()); // First component samples
    EXPECT_EQ(numPoints, samples[1].size()); // Second component samples

    // Check first component (circle) samples
    EXPECT_NEAR(1.0, samples[0][0].real(), 1e-10);
    EXPECT_NEAR(0.0, samples[0][0].imag(), 1e-10);

    // Check second component (ellipse) samples
    EXPECT_NEAR(2.0, samples[1][0].real(), 1e-10);
    EXPECT_NEAR(0.0, samples[1][0].imag(), 1e-10);
}

TEST_F(BoundaryTest, FindParameterizationDefaultComponent)
{
    Boundary boundary(m_circleComponent);

    Complex point(0.0, 1.0); // corresponds to t = π/2 on unit circle
    double param = boundary.findParameterization(point);
    EXPECT_NEAR(M_PI / 2, param, 1e-10);
}

TEST_F(BoundaryTest, FindParameterizationSpecificComponent)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    // Test on first component (circle)
    Complex point1(0.0, 1.0);
    double param1 = boundary.findParameterization(point1, 0);
    EXPECT_NEAR(M_PI / 2, param1, 1e-10);

    // Test on second component (ellipse with a=2.0, b=1.5)
    Complex point2(0.0, 1.5);
    double param2 = boundary.findParameterization(point2, 1);
    EXPECT_NEAR(M_PI / 2, param2, 1e-10);
}

TEST_F(BoundaryTest, GetComponent)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    const BoundaryComponent& comp1 = boundary.getComponent(0);
    const BoundaryComponent& comp2 = boundary.getComponent(1);

    // Verify we get the correct components by testing evaluation
    Complex pt1 = comp1.evaluate(0.0);
    EXPECT_NEAR(1.0, pt1.real(), 1e-10);

    Complex pt2 = comp2.evaluate(0.0);
    EXPECT_NEAR(2.0, pt2.real(), 1e-10);
}

TEST_F(BoundaryTest, GetComponents)
{
    std::vector<std::shared_ptr<BoundaryComponent>> components = {m_circleComponent, m_ellipseComponent};
    Boundary boundary(components);

    const auto& retrievedComponents = boundary.getComponents();
    EXPECT_EQ(2, retrievedComponents.size());
    EXPECT_EQ(m_circleComponent, retrievedComponents[0]);
    EXPECT_EQ(m_ellipseComponent, retrievedComponents[1]);
}

TEST_F(BoundaryTest, OutOfBoundsAccess)
{
    Boundary boundary(m_circleComponent);

    // These should throw exceptions for out-of-bounds access
    EXPECT_THROW(boundary.evaluate(0.0, 1), std::out_of_range);
    EXPECT_THROW(boundary.evaluateDerivative(0.0, 1), std::out_of_range);
    EXPECT_THROW(boundary.findParameterization(Complex(1.0, 0.0), 1), std::out_of_range);
    EXPECT_THROW(boundary.getComponent(1), std::out_of_range);
}
