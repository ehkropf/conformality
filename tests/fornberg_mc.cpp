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