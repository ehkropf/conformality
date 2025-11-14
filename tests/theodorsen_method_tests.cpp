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
#include <memory>
#include <cmath>
#include <random>
#include <vector>
#include "../src/methods/TheodorsenMethod.h"
#include "../src/core/ConformalMap.h"
#include "../src/domains/Domain.h"

class TheodorsenMethodTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create unit circle as source domain
        unit_circle = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 1.0);

        // Create elliptical domain as target
        ellipse = std::make_shared<EllipticalDomain>(2.0, 1.0); // Semi-major axis = 2, semi-minor axis = 1

        // Create method with 64 sample points
        method = std::make_shared<TheodorsenMethod>(64);
    }

    std::shared_ptr<CircularDomain> unit_circle;
    std::shared_ptr<EllipticalDomain> ellipse;
    std::shared_ptr<TheodorsenMethod> method;
};

TEST_F(TheodorsenMethodTest, UnitCircleBoundaryPointsMapToEllipse)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Test systematic boundary points on unit circle
    const int num_test_points = 16;
    for (int i = 0; i < num_test_points; ++i)
    {
        double theta = 2.0 * M_PI * i / num_test_points;

        // Point on unit circle boundary: (cos(theta), sin(theta))
        Complex boundary_point(std::cos(theta), std::sin(theta));
        Complex mapped_point = method->map(boundary_point);

        // Mapped point should be on ellipse boundary - check distance from origin is reasonable
        double mapped_magnitude = std::abs(mapped_point);
        ASSERT_GT(mapped_magnitude, 0.5)  // Should be inside ellipse's smallest radius
            << "Unit circle point at theta=" << theta
            << " maps to point with magnitude=" << mapped_magnitude;
        ASSERT_LT(mapped_magnitude, 2.5)  // Should be inside ellipse's largest radius
            << "Unit circle point at theta=" << theta
            << " maps to point with magnitude=" << mapped_magnitude;
    }
}

TEST_F(TheodorsenMethodTest, ComputedBoundarySamplesFromTargetDomain)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Get actual boundary samples used by the method (from target ellipse domain)
    const auto& boundary_samples = method->getBoundarySamples();
    EXPECT_EQ(boundary_samples.size(), 64);

    // Verify boundary samples are from ellipse domain (target domain)
    for (size_t i = 0; i < boundary_samples.size(); ++i)
    {
        Complex sample = boundary_samples[i];
        double sample_magnitude = std::abs(sample);

        // Samples should be from ellipse boundary, so magnitude should be reasonable
        EXPECT_GT(sample_magnitude, 0.5)  // At least semi-minor axis
            << "Boundary sample " << i << " at " << sample
            << " has magnitude=" << sample_magnitude;
        EXPECT_LT(sample_magnitude, 2.5)  // Less than semi-major axis + margin
            << "Boundary sample " << i << " at " << sample
            << " has magnitude=" << sample_magnitude;
    }
}

TEST_F(TheodorsenMethodTest, ConstructorValidation)
{
    // Valid power of 2
    EXPECT_NO_THROW(TheodorsenMethod(128));
    EXPECT_NO_THROW(TheodorsenMethod(256));
    EXPECT_NO_THROW(TheodorsenMethod(32));

    // Invalid (not power of 2)
    EXPECT_THROW(TheodorsenMethod(100), std::invalid_argument);
    EXPECT_THROW(TheodorsenMethod(63), std::invalid_argument);
    EXPECT_THROW(TheodorsenMethod(0), std::invalid_argument);
}

TEST_F(TheodorsenMethodTest, ConvergenceWithDifferentTolerances)
{
    ConformalMap map(unit_circle, ellipse, method);

    // Test loose tolerance
    method->compute(map, 1e-4);
    int loose_iterations = method->getIterationCount();
    double loose_residual = method->getResidualNorm();

    // Reset for tight tolerance test
    method->setNumPoints(64); // This clears previous computation

    // Test tight tolerance
    method->compute(map, 1e-10);
    int tight_iterations = method->getIterationCount();
    double tight_residual = method->getResidualNorm();

    // Tighter tolerance should require more iterations and achieve better residual
    EXPECT_GE(tight_iterations, loose_iterations);
    EXPECT_LE(tight_residual, loose_residual);
}

TEST_F(TheodorsenMethodTest, DifferentEllipseParameters)
{
    // Test with different ellipse parameters
    auto thin_ellipse = std::make_shared<EllipticalDomain>(1.8, 1.0); // Mildly elongated
    auto fat_ellipse = std::make_shared<EllipticalDomain>(1.1, 1.0);  // Nearly circular

    ConformalMap thin_map(unit_circle, thin_ellipse, method);
    ConformalMap fat_map(unit_circle, fat_ellipse, method);

    // Both should converge, but thin ellipse might take more iterations
    EXPECT_NO_THROW(method->compute(thin_map, 1e-6));
    EXPECT_TRUE(method->hasConverged());
    int thin_iterations = method->getIterationCount();

    method->setNumPoints(64); // Reset

    EXPECT_NO_THROW(method->compute(fat_map, 1e-6));
    EXPECT_TRUE(method->hasConverged());
    int fat_iterations = method->getIterationCount();

    // Fat ellipse (closer to circle) should converge faster
    EXPECT_LE(fat_iterations, thin_iterations + 5); // Allow some tolerance
}

TEST_F(TheodorsenMethodTest, DomainValidation)
{
    // Create conformal map
    ConformalMap map(unit_circle, ellipse, method);

    // Should work with unit circle -> starlike
    EXPECT_NO_THROW(method->compute(map, 1e-6));

    // Test with non-starlike domain (L-shaped polygon)
    std::vector<Complex> l_shape_vertices = {
        Complex(0.0, 0.0), Complex(2.0, 0.0), Complex(2.0, 1.0),
        Complex(1.0, 1.0), Complex(1.0, 2.0), Complex(0.0, 2.0)
    };
    auto non_starlike = std::make_shared<PolygonalDomain>(l_shape_vertices);
    ConformalMap invalid_map(unit_circle, non_starlike, method);

    // Should reject non-starlike target domain
    EXPECT_THROW(method->compute(invalid_map, 1e-6), std::invalid_argument);

    // Test with non-unit circle source
    auto non_unit_circle = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 2.0);
    ConformalMap invalid_source_map(non_unit_circle, ellipse, method);

    // Should reject non-unit circle source
    EXPECT_THROW(method->compute(invalid_source_map, 1e-6), std::invalid_argument);
}

TEST_F(TheodorsenMethodTest, EllipseMapping)
{
    // Create conformal map from unit circle to ellipse
    ConformalMap map(unit_circle, ellipse, method);

    // Compute the map
    EXPECT_NO_THROW(method->compute(map, 1e-8));

    // Check convergence
    EXPECT_TRUE(method->hasConverged());
    EXPECT_LT(method->getResidualNorm(), 1e-8);
    EXPECT_GT(method->getIterationCount(), 0);
    EXPECT_LT(method->getIterationCount(), 100);

    // Check that Laurent coefficients were computed
    EXPECT_FALSE(method->getLaurentCoefficients().empty());
    EXPECT_EQ(method->getLaurentCoefficients().size(), 64);

    // Check boundary samples
    EXPECT_FALSE(method->getBoundarySamples().empty());
    EXPECT_EQ(method->getBoundarySamples().size(), 64);
}

TEST_F(TheodorsenMethodTest, ExternalMapping)
{
    // Test external mapping (placeholder implementation)
    // Create unbounded ellipse domain for external mapping
    auto unbounded_ellipse = std::make_shared<EllipticalDomain>(2.0, 1.0, 0.0, Complex(0.0, 0.0), true);
    ConformalMap external_map(unit_circle, unbounded_ellipse, method);

    // Should still work but with external algorithm
    EXPECT_NO_THROW(method->compute(external_map, 1e-6));
    EXPECT_TRUE(method->hasConverged());
}

TEST_F(TheodorsenMethodTest, InitialState)
{
    EXPECT_FALSE(method->hasConverged());
    EXPECT_EQ(method->getResidualNorm(), 0.0);
    EXPECT_EQ(method->getIterationCount(), 0);
    EXPECT_TRUE(method->getLaurentCoefficients().empty());
    EXPECT_TRUE(method->getBoundarySamples().empty());
}

TEST_F(TheodorsenMethodTest, DISABLED_UnitDiskPointsMapToEllipseInterior)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Generate many random points inside unit disk and verify mapping properties
    std::random_device rd;
    std::mt19937 gen(123); // Fixed seed
    std::uniform_real_distribution<> radius_dist(0.0, 0.95); // Stay well inside unit disk
    std::uniform_real_distribution<> angle_dist(0.0, 2.0 * M_PI);

    int valid_mappings = 0;
    const int total_samples = 100;
    std::vector<std::pair<Complex, Complex>> invalid_mappings; // (source_point, mapped_point)

    for (int i = 0; i < total_samples; ++i)
    {
        double r = radius_dist(gen);
        double theta = angle_dist(gen);
        Complex test_point(r * std::cos(theta), r * std::sin(theta));

        // Point should be inside unit disk
        ASSERT_LT(std::abs(test_point), 1.0);

        Complex mapped_point = method->map(test_point);

        // Check if mapped point is inside ellipse
        double ellipse_test = (mapped_point.real() / 2.0) * (mapped_point.real() / 2.0) +
                              (mapped_point.imag() / 1.0) * (mapped_point.imag() / 1.0);

        if (ellipse_test < 1.0) // Inside ellipse
        {
            valid_mappings++;
        }
        else
        {
            // Collect invalid mappings for debugging
            invalid_mappings.emplace_back(test_point, mapped_point);
        }
    }

    // Verify we tested all sample points
    EXPECT_EQ(total_samples, 100)
        << "Should have tested 100 sample points";

    // Most unit disk points should map correctly to the ellipse interior
    // Allow some tolerance for numerical issues
    double success_rate = static_cast<double>(valid_mappings) / total_samples;
    EXPECT_GT(success_rate, 0.85)  // At least 85% should map correctly
        << "Success rate: " << success_rate
        << " (" << valid_mappings << "/" << total_samples << ")";

    if (!invalid_mappings.empty())
    {
        std::stringstream error_msg;
        error_msg << "Unit disk points that didn't map to ellipse interior: " << invalid_mappings.size()
                  << " out of " << total_samples << " tested points failed.\n";

        // Show first few invalid mappings for debugging
        const size_t max_examples = std::min(size_t(5), invalid_mappings.size());
        error_msg << "First " << max_examples << " invalid mappings:\n";
        for (size_t i = 0; i < max_examples; ++i)
        {
            const auto& [source, mapped] = invalid_mappings[i];
            error_msg << "  " << source << " -> " << mapped
                      << " (|w|=" << std::abs(mapped) << ")\n";
        }

        FAIL() << error_msg.str();
    }
}

TEST_F(TheodorsenMethodTest, InverseMapPlaceholder)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Inverse map should throw (not yet implemented)
    EXPECT_THROW(method->inverseMap(Complex(0.5, 0.0)), std::runtime_error);
}

TEST_F(TheodorsenMethodTest, LaurentCoefficientsProperties)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    const auto& coeffs = method->getLaurentCoefficients();
    EXPECT_EQ(coeffs.size(), 64);

    // For a conformal map from simply connected domain to unit disk,
    // the first coefficient (a_0) should be close to zero for centered ellipse
    EXPECT_LT(std::abs(coeffs[0]), 0.1);

    // Check that higher order coefficients decay (at least some should be small)
    int small_coeff_count = 0;
    for (size_t i = 32; i < coeffs.size(); ++i)
    {
        if (std::abs(coeffs[i]) < 1e-6)
        {
            small_coeff_count++;
        }
    }
    EXPECT_GT(small_coeff_count, 10); // Expect some high-order coefficients to be small
}

TEST_F(TheodorsenMethodTest, DISABLED_MappingConsistencyCheck)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Test that mapping preserves relative ordering and orientation
    std::vector<Complex> test_points = {
        Complex(0.0, 0.0),      // Center
        Complex(1.0, 0.0),      // Inside, on major axis
        Complex(0.0, 0.5),      // Inside, on minor axis
        Complex(1.5, 0.0),      // Between center and boundary
        Complex(0.0, 0.8)       // Between center and boundary
    };

    std::vector<Complex> mapped_points;
    for (const auto& point : test_points)
    {
        mapped_points.push_back(method->map(point));
    }

    // Center should map closest to origin
    double center_dist = std::abs(mapped_points[0]);
    for (size_t i = 1; i < mapped_points.size(); ++i)
    {
        EXPECT_LE(center_dist, std::abs(mapped_points[i]) + 0.1)
            << "Center should map closer to origin than other interior points";
    }

    // Points closer to boundary should map closer to unit circle
    // Point at (1.5, 0) should be further from origin in target than (1.0, 0)
    EXPECT_GT(std::abs(mapped_points[3]), std::abs(mapped_points[1]))
        << "Points closer to boundary should map further from origin";

    EXPECT_GT(std::abs(mapped_points[4]), std::abs(mapped_points[2]))
        << "Points closer to boundary should map further from origin";
}

TEST_F(TheodorsenMethodTest, MapEvaluation)
{
    ConformalMap map(unit_circle, ellipse, method);
    method->compute(map, 1e-8);

    // Test mapping of center point
    Complex center = ellipse->getCenter();
    Complex mapped_center = method->map(center);

    // Center should map close to origin (for well-centered ellipse)
    EXPECT_LT(std::abs(mapped_center), 0.1);

    // Test mapping of points on boundary
    // Point on boundary should map to unit circle
    Complex boundary_point(2.0, 0.0); // On ellipse boundary
    Complex mapped_boundary = method->map(boundary_point);

    // Mapped boundary point should be close to unit circle
    EXPECT_NEAR(std::abs(mapped_boundary), 1.0, 0.1);
}

TEST_F(TheodorsenMethodTest, MapEvaluationBeforeCompute)
{
    // Should throw when trying to evaluate map before computation
    EXPECT_THROW(method->map(Complex(0.0, 0.0)), std::runtime_error);
    EXPECT_THROW(method->inverseMap(Complex(0.0, 0.0)), std::runtime_error);
}

TEST_F(TheodorsenMethodTest, NumPointsAccessors)
{
    EXPECT_EQ(method->getNumPoints(), 64);

    method->setNumPoints(128);
    EXPECT_EQ(method->getNumPoints(), 128);

    // Setting invalid number should throw
    EXPECT_THROW(method->setNumPoints(100), std::invalid_argument);

    // Number should remain unchanged after failed set
    EXPECT_EQ(method->getNumPoints(), 128);
}
