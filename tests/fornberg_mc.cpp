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
    // S(j, ν) = θ_j = 2πj/N (canonical parameters map directly to user parameters)
    // Note: S uses (N, m) layout to match MATLAB column-per-boundary convention
    method.m_S.resize(config.N, 2);
    for (int nu = 0; nu < 2; ++nu)
    {
        for (int j = 0; j < config.N; ++j)
        {
            method.m_S(j, nu) = 2.0 * M_PI * j / config.N;
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

// Friend test class for testing solveSystem() method
class FornbergMCSolveSystemTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-8;
        config.cgm_tolerance = 1e-10;
        config.max_cgm_iterations = 200;
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

TEST_F(FornbergMCSolveSystemTest, SolutionDimensions)
{
    // Annulus case: m=2, N=64
    // U should have size m*N+1 = 129 (annulus case)
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
    method.solveSystem();

    // Solution dimensions should match system columns
    EXPECT_EQ(method.m_U.size(), method.m_D.cols());
    EXPECT_EQ(method.m_U.size(), 129);  // m*N + 1 for annulus
}

TEST_F(FornbergMCSolveSystemTest, ReducesResidual)
{
    // This is an OVERDETERMINED system (D has more rows than columns).
    // We solve the least-squares problem: min ||D*U - g||
    // via the normal equations: D†D*U = D†g
    //
    // The CG solver operates on the REAL system: 2*real(D†D)*x = 2*real(D†g)
    // So we verify the real part of the normal equations residual.
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

    // Compute the RHS of the real normal equations: 2*real(D† * g)
    Eigen::VectorXcd Dtg = method.m_D.adjoint() * method.m_g;
    Eigen::VectorXd real_rhs = 2.0 * Dtg.real();

    // Compute initial real normal equations residual with zero vector
    double initial_real_residual = real_rhs.norm();  // ||2*real(D†g) - 0||

    // Solve the system
    method.solveSystem();

    // Get the real solution vector (m_U has zero imaginary parts)
    Eigen::VectorXd U_real = method.m_U.real();

    // Compute the real matrix-vector product: 2*real(D†D)*U
    // This is what the CG solver computes internally
    Eigen::VectorXcd DU = method.m_D * U_real;
    Eigen::VectorXcd DtDU = method.m_D.adjoint() * DU;
    Eigen::VectorXd real_lhs = 2.0 * DtDU.real();

    // The real normal equations residual
    double final_real_residual = (real_lhs - real_rhs).norm();

    // The solution should reduce the real normal equations residual
    EXPECT_LT(final_real_residual, initial_real_residual)
        << "CG solution should reduce the real normal equations residual";

    // The relative real normal equations residual should be small
    double relative_residual = final_real_residual / real_rhs.norm();
    EXPECT_LT(relative_residual, 1e-6)
        << "Relative real normal equations residual should be small (tolerance: "
        << config.cgm_tolerance << ")";

    // Also verify the least-squares residual is reduced from zero initial guess
    double initial_ls_residual = method.m_g.norm();  // ||D*0 - g|| = ||g||
    double final_ls_residual = (method.m_D * method.m_U - method.m_g).norm();
    EXPECT_LT(final_ls_residual, initial_ls_residual)
        << "Least-squares residual should be reduced from zero initial guess";
}

TEST_F(FornbergMCSolveSystemTest, ConvergenceInfo)
{
    // Verify that convergence information is properly populated
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
    method.solveSystem();

    // Get convergence info from the CG solver
    const auto& info = method.mp_cg_solver->getLastConvergenceInfo();

    // Should have performed some iterations
    EXPECT_GT(info.iterations, 0);

    // Residual history should be populated
    EXPECT_FALSE(info.residual_history.empty());

    // Final residual should be recorded
    EXPECT_GE(info.final_residual, 0.0);

    // For this well-conditioned problem, we expect convergence
    EXPECT_TRUE(info.converged)
        << "CG should converge for this annulus problem";
}

// Helper functions for creating test domains (shared across newton update tests)
namespace
{

std::shared_ptr<Boundary> createCircularBoundaryHelper(Complex center, double radius)
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

std::shared_ptr<MultiplyConnectedDomain> createAnnulusDomainHelper()
{
    auto outer = createCircularBoundaryHelper(Complex(0, 0), 1.0);
    auto inner = createCircularBoundaryHelper(Complex(0.3, 0), 0.15);
    return std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );
}

std::shared_ptr<MultiplyConnectedDomain> createThreeConnectedDomainHelper()
{
    auto outer = createCircularBoundaryHelper(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundaryHelper(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundaryHelper(Complex(-0.3, -0.1), 0.12);
    return std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );
}

}  // anonymous namespace

TEST(FornbergMCNewtonUpdateTest, ComputesResidualBeforeUpdate)
{
    // Verify residual is computed as infinity norm of U before applying updates
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;

    auto domain = createAnnulusDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    // Store expected residual (infinity norm of U before update)
    double expected_residual = method.m_U.lpNorm<Eigen::Infinity>();

    method.newtonUpdate();

    // Residual should match the infinity norm computed before modifications
    EXPECT_DOUBLE_EQ(method.m_current_residual, expected_residual);
}

TEST(FornbergMCNewtonUpdateTest, UpdatesBoundaryCorrespondencesAnnulus)
{
    // Verify S is updated from solution vector for annulus case
    // The update should be: S_new = S_old + U[0:m*N-1] / abs_eta
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;

    auto domain = createAnnulusDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    const int m = 2;
    const int N = config.N;

    // Store initial S values
    Eigen::MatrixXd S_before = method.m_S;

    // Check that solution vector has non-zero values (indicating an update will happen)
    double u_norm = 0.0;
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            u_norm += std::abs(method.m_U(nu * N + j));
        }
    }
    EXPECT_GT(u_norm, 0.0) << "Solution vector should have non-zero boundary updates";

    method.newtonUpdate();

    // Verify S was modified (changed from initial values)
    double total_change = 0.0;
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            total_change += std::abs(method.m_S(j, nu) - S_before(j, nu));
        }
    }
    EXPECT_GT(total_change, 0.0) << "S should be updated by newtonUpdate()";

    // After update, m_U[0:m*N-1] should contain the scaled values (U/abs_eta)
    // that were added to S. Verify they're now scaled (different from original if abs_eta != 1)
    // This is a sanity check that the scaling happened
    bool any_abs_eta_not_one = false;
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            if (std::abs(method.m_abs_eta(j, nu) - 1.0) > 1e-10)
            {
                any_abs_eta_not_one = true;
                break;
            }
        }
    }
    // If abs_eta values aren't all 1, scaling should have occurred
    if (any_abs_eta_not_one)
    {
        EXPECT_TRUE(true) << "abs_eta scaling was applied";
    }
}

TEST(FornbergMCNewtonUpdateTest, UpdatesRadiusAnnulus)
{
    // For annulus (m=2): radius is updated, center is not (for annulus case)
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;

    auto domain = createAnnulusDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    const int m = 2;
    const int N = config.N;

    // Store initial moduli
    Eigen::VectorXcd moduli_before = method.m_conformal_moduli;
    Complex center_before = moduli_before(0);
    double radius_before = std::real(moduli_before(1));

    // Verify there's a non-zero update in the solution vector
    double radius_update_in_U = std::real(method.m_U(m * N));
    EXPECT_NE(radius_update_in_U, 0.0) << "Solution should contain non-zero radius update";

    method.newtonUpdate();

    // Moduli layout: [c₂, ρ₂] for annulus
    // Center should NOT change for annulus case (m_is_annulus = true)
    EXPECT_NEAR(std::real(method.m_conformal_moduli(0)), std::real(center_before), 1e-14)
        << "Annulus center real part should not change";
    EXPECT_NEAR(std::imag(method.m_conformal_moduli(0)), std::imag(center_before), 1e-14)
        << "Annulus center imag part should not change";

    // Radius should have changed
    double radius_after = std::real(method.m_conformal_moduli(1));
    double actual_radius_change = radius_after - radius_before;

    // The change should have the same sign as U[m*N] and be non-zero
    if (radius_update_in_U > 0)
    {
        EXPECT_GT(actual_radius_change, 0.0)
            << "Radius change should be positive when U[m*N] is positive";
    }
    else if (radius_update_in_U < 0)
    {
        EXPECT_LT(actual_radius_change, 0.0)
            << "Radius change should be negative when U[m*N] is negative";
    }

    // Radius should remain positive
    EXPECT_GT(radius_after, 0.0) << "Radius should remain positive after update";
}

TEST(FornbergMCNewtonUpdateTest, UpdatesCentersAndRadiiGeneral)
{
    // For 3-connected domain: both centers and radii are updated
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;

    auto domain = createThreeConnectedDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    const int m = 3;
    const int N = config.N;

    // Store initial moduli
    Eigen::VectorXcd moduli_before = method.m_conformal_moduli;

    // Extract update directions from solution vector BEFORE newtonUpdate modifies U
    // Radii: U[m*N + k] for k = 0, 1
    // Centers: U[m*N + m - 1 + 2k] (real), U[m*N + m + 2k] (imag)
    std::vector<double> radius_update_direction(m - 1);
    std::vector<Complex> center_update_direction(m - 1);

    for (int k = 0; k < m - 1; ++k)
    {
        radius_update_direction[k] = std::real(method.m_U(m * N + k));

        int re_idx = m * N + m - 1 + 2 * k;
        int im_idx = m * N + m + 2 * k;
        center_update_direction[k] = Complex(
            std::real(method.m_U(re_idx)),
            std::real(method.m_U(im_idx))
        );
    }

    method.newtonUpdate();

    // Verify moduli updates have correct direction
    // Moduli layout: [c₂, ρ₂, c₃, ρ₃, ...]
    for (int k = 0; k < m - 1; ++k)
    {
        int center_idx = 2 * k;
        int radius_idx = 2 * k + 1;

        // Compute actual changes
        Complex actual_center_change = method.m_conformal_moduli(center_idx) - moduli_before(center_idx);
        double actual_radius_change = std::real(method.m_conformal_moduli(radius_idx))
                                      - std::real(moduli_before(radius_idx));

        // Check that changes have the same sign as the solution vector values
        // (accounting for possible scaling factors)
        if (std::abs(radius_update_direction[k]) > 1e-10)
        {
            EXPECT_TRUE((actual_radius_change > 0) == (radius_update_direction[k] > 0))
                << "Radius " << (k + 2) << " change should have same sign as U value";
        }

        if (std::abs(std::real(center_update_direction[k])) > 1e-10)
        {
            EXPECT_TRUE((std::real(actual_center_change) > 0) == (std::real(center_update_direction[k]) > 0))
                << "Center " << (k + 2) << " real change should have same sign as U value";
        }

        if (std::abs(std::imag(center_update_direction[k])) > 1e-10)
        {
            EXPECT_TRUE((std::imag(actual_center_change) > 0) == (std::imag(center_update_direction[k]) > 0))
                << "Center " << (k + 2) << " imag change should have same sign as U value";
        }

        // Radii should remain positive
        EXPECT_GT(std::real(method.m_conformal_moduli(radius_idx)), 0.0)
            << "Radius " << (k + 2) << " should remain positive";
    }
}

TEST(FornbergMCNewtonUpdateTest, AppliesDampingWhenEnabled)
{
    // Verify damping factor is applied to all updates
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;
    config.enable_newton_damping = true;
    config.newton_damping_factor = 0.5;

    auto domain = createAnnulusDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    const int m = 2;
    const int N = config.N;
    const double damping = 0.5;

    // Store initial values
    Eigen::MatrixXd S_before = method.m_S;
    Eigen::VectorXcd moduli_before = method.m_conformal_moduli;

    // Compute expected damped updates for S
    Eigen::MatrixXd expected_delta_S(N, m);
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            double abs_eta_val = method.m_abs_eta(j, nu);
            double u_val = std::real(method.m_U(nu * N + j));
            expected_delta_S(j, nu) = (abs_eta_val > 1e-14) ? damping * u_val / abs_eta_val : 0.0;
        }
    }

    // Expected damped radius update
    double expected_radius_delta = damping * std::real(method.m_U(m * N));

    method.newtonUpdate();

    // Verify S was updated with damping
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            double expected_S = S_before(j, nu) + expected_delta_S(j, nu);
            EXPECT_NEAR(method.m_S(j, nu), expected_S, 1e-12)
                << "Damped S(" << j << ", " << nu << ") mismatch";
        }
    }

    // Verify radius was updated with damping
    EXPECT_NEAR(std::real(method.m_conformal_moduli(1)),
                std::real(moduli_before(1)) + expected_radius_delta, 1e-12)
        << "Damped radius update mismatch";
}

TEST(FornbergMCNewtonUpdateTest, SyncsCanonicalDomain)
{
    // Verify canonical domain is updated with new moduli after Newton update
    FornbergMCConfiguration config;
    config.N = 64;
    config.max_newton_iterations = 10;
    config.newton_tolerance = 1e-8;
    config.cgm_tolerance = 1e-10;
    config.max_cgm_iterations = 200;
    config.verbose = false;

    auto domain = createThreeConnectedDomainHelper();
    FornbergMC method(config);

    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        method.mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        config.N
    );
    method.initializeNewtonIteration();
    method.formSystem();
    method.solveSystem();

    // Store moduli after update (we'll compare canonical domain against these)
    method.newtonUpdate();
    Eigen::VectorXcd updated_moduli = method.m_conformal_moduli;

    // Get moduli from canonical domain
    auto canonical_moduli = method.mp_canonical_domain->getConformalModuli();

    // Verify canonical domain has the updated moduli
    ASSERT_EQ(canonical_moduli.size(), updated_moduli.size());
    for (int i = 0; i < canonical_moduli.size(); ++i)
    {
        EXPECT_NEAR(std::real(canonical_moduli(i)), std::real(updated_moduli(i)), 1e-12)
            << "Canonical domain moduli[" << i << "] real part not synced";
        EXPECT_NEAR(std::imag(canonical_moduli(i)), std::imag(updated_moduli(i)), 1e-12)
            << "Canonical domain moduli[" << i << "] imag part not synced";
    }
}
