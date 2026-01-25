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
#include "../src/methods/FornbergMC.h"
#include "../src/methods/FornbergMCConfiguration.h"
#include "../src/methods/PMatrixBuilder.h"
#include "../src/numerics/CGSolver.h"
#include "../src/domains/FornbergCanonicalDomain.h"
#include "../src/domains/Domain.h"
#include "../src/core/ConformalMap.h"

class FornbergMCTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create test configuration
        config.N = 64; // Smaller N for faster tests
        config.max_newton_iterations = 10;
        config.max_cgm_iterations = 100;
        config.newton_tolerance = 1e-8;
        config.cgm_tolerance = 1e-8;
        config.verbose = false; // Quiet tests
    }

    FornbergMCConfiguration config;
};

// Test FornbergMCConfiguration validation
TEST_F(FornbergMCTest, ConfigurationValidation)
{
    // Valid configuration should pass
    EXPECT_NO_THROW(config.validate());

    // Newton parameter validation
    config.newton_tolerance = -1.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.newton_tolerance = 1e-17; // Too small
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.newton_tolerance = 1e-8; // Reset

    // Newton damping factor validation
    config.newton_damping_factor = 0.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.newton_damping_factor = 1.5;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.newton_damping_factor = 0.8; // Reset

    // CG parameter validation
    config.cgm_tolerance = -1.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.cgm_tolerance = 1e-17; // Too small
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.cgm_tolerance = 1e-8; // Reset

    config.cgm_restart_threshold = 0.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.cgm_restart_threshold = 1.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.cgm_restart_threshold = 0.1; // Reset

    // Tolerance consistency
    config.cgm_tolerance = 1e-6;
    config.newton_tolerance = 1e-8; // CG more relaxed than Newton
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.cgm_tolerance = 1e-8; // Reset
    config.newton_tolerance = 1e-8; // Reset

    // Discretization validation
    config.N = 63; // Not power of 2
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.N = 64; // Reset

    config.max_N = 32; // Less than N
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.max_N = 128; // Reset

    // Adaptive N with non-power-of-2 max_N
    config.adaptive_N = true;
    config.max_N = 100; // Not power of 2
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.max_N = 128; // Reset
    config.adaptive_N = false; // Reset

    // Redistribution validation
    config.redistribution_threshold = 0.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.redistribution_threshold = 1.5;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.redistribution_threshold = 0.1; // Reset

    // Annulus threshold validation
    config.annulus_aspect_ratio_threshold = 1.0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.annulus_aspect_ratio_threshold = 10.0; // Reset

    // Thread count validation
    config.parallel_threads = 0;
    EXPECT_THROW(config.validate(), std::invalid_argument);
    config.parallel_threads = 1; // Reset
}

TEST_F(FornbergMCTest, ConfigurationPresets)
{
    FornbergMCConfiguration high_precision_config;
    high_precision_config.setHighPrecision();
    EXPECT_EQ(high_precision_config.newton_tolerance, 1e-15);
    EXPECT_EQ(high_precision_config.cgm_tolerance, 1e-15);
    EXPECT_TRUE(high_precision_config.monitor_eigenvalues);

    FornbergMCConfiguration fast_config;
    fast_config.setFastComputation();
    EXPECT_EQ(fast_config.newton_tolerance, 1e-8);
    EXPECT_EQ(fast_config.max_newton_iterations, 20);
    EXPECT_FALSE(fast_config.verbose);

    FornbergMCConfiguration annulus_config;
    annulus_config.setAnnulusOptimized();
    EXPECT_TRUE(annulus_config.auto_detect_annulus);
    EXPECT_TRUE(annulus_config.enable_best_iterate);
}

// Test FornbergMC construction
TEST_F(FornbergMCTest, Construction)
{
    // Default construction should work
    EXPECT_NO_THROW(FornbergMC method);

    // Construction with valid config should work
    EXPECT_NO_THROW(FornbergMC method(config));

    // Construction with invalid config should fail
    config.N = 63; // Not power of 2
    EXPECT_THROW(FornbergMC method(config), std::invalid_argument);
}

TEST_F(FornbergMCTest, ConfigurationAccess)
{
    FornbergMC method(config);

    // Should be able to get configuration
    const auto& retrieved_config = method.getConfiguration();
    EXPECT_EQ(retrieved_config.N, config.N);
    EXPECT_EQ(retrieved_config.newton_tolerance, config.newton_tolerance);

    // Should be able to update configuration
    FornbergMCConfiguration new_config = config;
    new_config.N = 128;
    EXPECT_NO_THROW(method.setConfiguration(new_config));
    EXPECT_EQ(method.getConfiguration().N, 128);
}

// Test domain validation - disabled because validation methods are protected
TEST_F(FornbergMCTest, DISABLED_DomainValidation)
{
    // This test is disabled because validateSourceDomain and validateTargetDomain
    // are protected methods. Domain validation will be tested through the compute() method
    // TODO: Create proper integration tests that exercise domain validation
}


// Test FornbergCanonicalDomain
TEST_F(FornbergMCTest, FornbergCanonicalDomain)
{
    // Test basic construction with hole parameters
    std::vector<Complex> centers = {Complex(0.3, 0.0)};
    std::vector<double> radii = {0.1};

    EXPECT_NO_THROW(FornbergCanonicalDomain canonical_domain(centers, radii, 64));

    FornbergCanonicalDomain canonical_domain(centers, radii, 64);

    // Test domain properties
    EXPECT_EQ(canonical_domain.getConnectivity(), 2);  // Unit disk + 1 hole
    EXPECT_TRUE(canonical_domain.isAnnulus());
    EXPECT_EQ(canonical_domain.getHoleCenters().size(), 1);
    EXPECT_EQ(canonical_domain.getHoleRadii().size(), 1);
    EXPECT_EQ(canonical_domain.getBoundaryPointCount(), 64);

    // Test conformal moduli
    auto moduli = canonical_domain.getConformalModuli();
    EXPECT_EQ(moduli.size(), 2);  // One pair (center, radius)

    // Test parameter updates
    std::vector<Complex> new_centers = {Complex(0.2, 0.1)};
    std::vector<double> new_radii = {0.08};
    EXPECT_NO_THROW(canonical_domain.updateHoleParameters(new_centers, new_radii));

    // Test validation
    EXPECT_TRUE(canonical_domain.isValidConfiguration());

    // Test invalid configurations
    std::vector<Complex> bad_centers = {Complex(0.0, 0.0)};
    std::vector<double> bad_radii = {1.1};  // Too large
    EXPECT_THROW(FornbergCanonicalDomain bad_domain(bad_centers, bad_radii, 64), std::invalid_argument);
}

// Integration test placeholder
TEST_F(FornbergMCTest, DISABLED_BasicComputeIntegration)
{
    // This test is disabled until we have proper multiply connected domains set up
    // TODO: Implement when domain creation utilities are available

    FornbergMC method(config);

    // Create source and target domains
    // auto source_domain = createCanonicalDomain();
    // auto target_domain = createTestMultiplyConnectedDomain();
    //
    // ConformalMap map(source_domain, target_domain);
    //
    // EXPECT_NO_THROW(method.compute(map));
}

// Convergence history test
TEST_F(FornbergMCTest, ConvergenceTracking)
{
    FornbergMC method(config);

    // Initially should not be converged
    EXPECT_FALSE(method.hasConverged());

    // Should have empty residual history
    EXPECT_TRUE(method.getResidualHistory().empty());

    // Current residual should be initialized to large value
    EXPECT_GT(method.getCurrentResidual(), 1e10);
}

// Test that Eigen integration works
TEST_F(FornbergMCTest, EigenIntegration)
{
    // Test that we can create and use Eigen matrices
    Eigen::MatrixXcd test_matrix = Eigen::MatrixXcd::Identity(3, 3);
    Eigen::VectorXcd test_vector = Eigen::VectorXcd::Ones(3);

    Eigen::VectorXcd result = test_matrix * test_vector;

    // Should get identity result
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(result[i].real(), 1.0, 1e-15);
        EXPECT_NEAR(result[i].imag(), 0.0, 1e-15);
    }
}

// Friend test class for testing formSystem() method
class FornbergMCFormSystemTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-8;
        config.cgm_tolerance = 1e-8;
        config.verbose = false;
    }

    // Create a circular boundary centered at the given point with given radius
    std::shared_ptr<Boundary> createCircularBoundary(Complex center, double radius)
    {
        auto component = std::make_shared<AnalyticBoundaryComponent>(
            [center, radius](double theta) {
                return center + radius * Complex(std::cos(theta), std::sin(theta));
            },
            [radius](double theta) {
                return radius * Complex(-std::sin(theta), std::cos(theta));
            }
        );
        return std::make_shared<Boundary>(component);
    }

    // Create an annulus domain (2-connected)
    std::shared_ptr<MultiplyConnectedDomain> createAnnulusDomain()
    {
        auto outer = createCircularBoundary(Complex(0, 0), 1.0);
        auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
        return std::make_shared<MultiplyConnectedDomain>(
            std::vector<std::shared_ptr<Boundary>>{outer, inner}
        );
    }

    // Create a 3-connected domain (unit disk with 2 holes)
    std::shared_ptr<MultiplyConnectedDomain> createThreeConnectedDomain()
    {
        auto outer = createCircularBoundary(Complex(0, 0), 1.0);
        auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
        auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
        return std::make_shared<MultiplyConnectedDomain>(
            std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
        );
    }

    FornbergMCConfiguration config;
};

TEST_F(FornbergMCFormSystemTest, DimensionsAnnulus)
{
    // Annulus (m=2, N=64)
    // Expected: D is (m*M x m*N+1) = (64 x 129)
    auto domain = createAnnulusDomain();

    FornbergMC method(config);

    // Set up internal state for testing
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;

    // Create canonical domain from user domain
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );

    method.initializeNewtonIteration();
    method.formSystem();

    EXPECT_EQ(method.getSystemMatrix().rows(), 64);
    EXPECT_EQ(method.getSystemMatrix().cols(), 129);
    EXPECT_EQ(method.getRHSVector().size(), 64);
}

TEST_F(FornbergMCFormSystemTest, DimensionsGeneral)
{
    // 3-connected domain (m=3, N=64)
    // Expected: D is (m*M+2 x m*N+3*(m-1)) = (98 x 198)
    auto domain = createThreeConnectedDomain();

    FornbergMC method(config);

    // Set up internal state for testing
    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;

    // Create canonical domain from user domain
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );

    method.initializeNewtonIteration();
    method.formSystem();

    EXPECT_EQ(method.getSystemMatrix().rows(), 98);
    EXPECT_EQ(method.getSystemMatrix().cols(), 198);
    EXPECT_EQ(method.getRHSVector().size(), 98);
}

TEST_F(FornbergMCFormSystemTest, NonZeroOutput)
{
    // Verify D and g are populated (not all zeros)
    auto domain = createAnnulusDomain();

    FornbergMC method(config);

    // Set up internal state for testing
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;

    // Create canonical domain from user domain
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );

    method.initializeNewtonIteration();
    method.formSystem();

    EXPECT_GT(method.getSystemMatrix().norm(), 0.0);
    EXPECT_GT(method.getRHSVector().norm(), 0.0);
}

// Friend test class for testing private computeFourierCoefficients() method
class FornbergMCFourierTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-8;
        config.verbose = false;
    }

    FornbergMCConfiguration config;
};

// Test computeFourierCoefficients() with a simple circular boundary
TEST_F(FornbergMCFourierTest, CircularBoundaryCoefficients)
{
    // Create a simple annulus: unit circle with a small circular hole
    // Outer boundary: circle at origin with radius 1.0
    // Inner boundary: circle at (0.3, 0) with radius 0.15

    auto outer_component = std::make_shared<AnalyticBoundaryComponent>(
        [](double theta) { return Complex(std::cos(theta), std::sin(theta)); },
        [](double theta) { return Complex(-std::sin(theta), std::cos(theta)); } // derivative
    );

    auto inner_component = std::make_shared<AnalyticBoundaryComponent>(
        [](double theta) {
            return Complex(0.3, 0.0) + 0.15 * Complex(std::cos(theta), std::sin(theta));
        },
        [](double theta) {
            return 0.15 * Complex(-std::sin(theta), std::cos(theta)); // derivative
        }
    );

    auto outer_boundary = std::make_shared<Boundary>(outer_component);
    auto inner_boundary = std::make_shared<Boundary>(inner_component);

    std::vector<std::shared_ptr<Boundary>> boundaries = {outer_boundary, inner_boundary};
    auto test_domain = std::make_shared<MultiplyConnectedDomain>(boundaries);

    // Create FornbergMC instance
    FornbergMC method(config);

    // Manually set up internal state to test computeFourierCoefficients()
    // This is normally done during compute(), but we're testing the method in isolation

    // Set connectivity
    method.m_connectivity = 2;

    // Set user domain
    method.mp_user_domain = test_domain;

    // Initialize S matrix with identity correspondence
    // S(ν, j) = θ_j = 2πj/N (canonical parameters map directly to user parameters)
    method.m_S.resize(2, config.N);
    for (int nu = 0; nu < 2; ++nu)
    {
        for (int j = 0; j < config.N; ++j)
        {
            method.m_S(nu, j) = 2.0 * M_PI * j / config.N;
        }
    }

    // Initialize coefficient matrix
    method.m_a.resize(config.N, 2);
    method.m_a.setZero();

    // Call the method we're testing
    method.computeFourierCoefficients();

    // Verify results
    const auto& coeffs = method.m_a;

    // Check dimensions
    EXPECT_EQ(coeffs.rows(), config.N);
    EXPECT_EQ(coeffs.cols(), 2);

    // Check that coefficients are not NaN or inf
    for (int nu = 0; nu < 2; ++nu)
    {
        for (int j = 0; j < config.N; ++j)
        {
            EXPECT_FALSE(std::isnan(coeffs(j, nu).real()))
                << "NaN found at coefficient (" << j << ", " << nu << ")";
            EXPECT_FALSE(std::isnan(coeffs(j, nu).imag()))
                << "NaN found at coefficient (" << j << ", " << nu << ")";
            EXPECT_FALSE(std::isinf(coeffs(j, nu).real()))
                << "Inf found at coefficient (" << j << ", " << nu << ")";
            EXPECT_FALSE(std::isinf(coeffs(j, nu).imag()))
                << "Inf found at coefficient (" << j << ", " << nu << ")";
        }
    }

    // For the outer boundary (unit circle at origin):
    // z(θ) = e^(iθ) = cos(θ) + i*sin(θ) has Fourier series with one dominant term
    // The FFT should have significant magnitude at index 1 (corresponding to e^(iθ))

    double max_coeff_magnitude = 0.0;
    int max_coeff_index = 0;
    for (int j = 0; j < config.N; ++j)
    {
        double mag = std::abs(coeffs(j, 0));
        if (mag > max_coeff_magnitude)
        {
            max_coeff_magnitude = mag;
            max_coeff_index = j;
        }
    }

    // The dominant coefficient should be at index 1 (positive frequency for e^(iθ))
    EXPECT_EQ(max_coeff_index, 1)
        << "Dominant Fourier coefficient should be at index 1 for e^(iθ) parameterization";

    // The magnitude should be approximately 1.0 (the radius)
    EXPECT_NEAR(max_coeff_magnitude, 1.0, 0.1)
        << "Dominant coefficient magnitude should match circle radius";

    // Sum of all coefficient magnitudes should be reasonable (not exploding)
    double total_magnitude = 0.0;
    for (int nu = 0; nu < 2; ++nu)
    {
        for (int j = 0; j < config.N; ++j)
        {
            total_magnitude += std::abs(coeffs(j, nu));
        }
    }
    EXPECT_LT(total_magnitude, 10.0)
        << "Total coefficient magnitude should be bounded for simple circular boundaries";
}
