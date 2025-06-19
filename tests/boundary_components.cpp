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
#include "../src/domains/BoundaryComponent.h"

TEST(BoundaryComponentTest, AnalyticEvaluation) {
    // Create an analytic boundary component for a circle
    auto circleFunc = [](double t) -> Complex {
        return Complex(std::cos(t), std::sin(t));
    };
    auto circleDerivFunc = [](double t) -> Complex {
        return Complex(-std::sin(t), std::cos(t));
    };

    AnalyticBoundaryComponent circle(circleFunc, circleDerivFunc);

    // Test evaluation at specific points
    Complex pt1 = circle.evaluate(0.0);
    EXPECT_NEAR(1.0, pt1.real(), 1e-10);
    EXPECT_NEAR(0.0, pt1.imag(), 1e-10);

    Complex pt2 = circle.evaluate(M_PI / 2);
    EXPECT_NEAR(0.0, pt2.real(), 1e-10);
    EXPECT_NEAR(1.0, pt2.imag(), 1e-10);

    // Test derivative
    Complex deriv1 = circle.evaluateDerivative(0.0);
    EXPECT_NEAR(0.0, deriv1.real(), 1e-10);
    EXPECT_NEAR(1.0, deriv1.imag(), 1e-10);
}

TEST(BoundaryComponentTest, Sampling) {
    // Create an analytic boundary component for an ellipse
    double a = 2.0, b = 1.0;
    auto ellipseFunc = [a, b](double t) -> Complex {
        return Complex(a * std::cos(t), b * std::sin(t));
    };
    auto ellipseDerivFunc = [a, b](double t) -> Complex {
        return Complex(-a * std::sin(t), b * std::cos(t));
    };

    AnalyticBoundaryComponent ellipse(ellipseFunc, ellipseDerivFunc);

    // Test sampling
    int numPoints = 4;
    auto samples = ellipse.sample(numPoints);

    EXPECT_EQ(numPoints, samples.size());

    // Check specific points (0, π/2, π, 3π/2)
    EXPECT_NEAR(a, samples[0].real(), 1e-10);
    EXPECT_NEAR(0.0, samples[0].imag(), 1e-10);

    EXPECT_NEAR(0.0, samples[1].real(), 1e-10);
    EXPECT_NEAR(b, samples[1].imag(), 1e-10);
}

TEST(BoundaryComponentTest, Parameterization) {
    // Create an analytic boundary component for a circle
    auto circleFunc = [](double t) -> Complex {
        return Complex(std::cos(t), std::sin(t));
    };
    auto circleDerivFunc = [](double t) -> Complex {
        return Complex(-std::sin(t), std::cos(t));
    };

    AnalyticBoundaryComponent circle(circleFunc, circleDerivFunc);

    // Find parameter for a point on the boundary
    Complex point(0.0, 1.0); // corresponds to t = π/2
    double param = circle.findParameterization(point);
    EXPECT_NEAR(M_PI / 2, param, 1e-10);
}
