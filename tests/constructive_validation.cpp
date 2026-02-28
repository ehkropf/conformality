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

#include "../src/core/StatusManager.h"

#include <cmath>

// These tests use FornbergMC's internal Newton iteration methods directly
// rather than the public compute() API. This enables fine-grained inspection
// of intermediate state (residuals, moduli, coefficients) at each stage,
// which is not exposed through the public API.

static std::shared_ptr<Boundary> createCircularBoundary(Complex center, double radius)
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

// ============================================================================
// Base Fixture
// ============================================================================

class ConstructiveValidationBase : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;

    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 20;
        config.newton_tolerance = 1e-14;
        config.cgm_tolerance = 1e-15;
        config.max_cgm_iterations = 20;
        config.verbose = false;
    }
};

// ============================================================================
// Identity Map Tests
// ============================================================================

class ConstructiveIdentityMap : public ConstructiveValidationBase {};

TEST_F(ConstructiveIdentityMap, IdentityM4ConvergesToTargetModuli)
{
    // th_gen_ex3 domain from Kropf thesis: 4-connected, all circular boundaries.
    // Since all user boundaries are circles, the conformal map from
    // the converged canonical domain to the user domain is the identity.
    // The algorithm must discover this by iterating the moduli from
    // an imprecise initial guess to match the actual boundary geometry.
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    // Initial guess: intentionally offset from actual centers to test convergence
    std::vector<Complex> hole_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
    std::vector<double> hole_radii = {0.25, 0.25, 0.25};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        hole_centers, hole_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Identity m=4 should converge within "
        << config.max_newton_iterations << " iterations";

    // Converged moduli must match target boundary centers and radii
    const std::vector<Complex> target_centers = {Complex(-0.5, 0.0), Complex(0.25, 0.43), Complex(0.25, -0.43)};
    const std::vector<double> target_radii = {0.25, 0.25, 0.25};

    const auto& actual_centers = method.mp_canonical_domain->getHoleCenters();
    const auto& actual_radii = method.mp_canonical_domain->getHoleRadii();

    for (size_t i = 0; i < target_centers.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual_centers[i] - target_centers[i]), 0.0, 1e-4)
            << "Center c[" << i << "]: got " << actual_centers[i]
            << ", expected " << target_centers[i];
    }
    for (size_t i = 0; i < target_radii.size(); ++i)
    {
        EXPECT_NEAR(actual_radii[i], target_radii[i], 1e-4)
            << "Radius rho[" << i << "]: got " << actual_radii[i]
            << ", expected " << target_radii[i];
    }
}

TEST_F(ConstructiveIdentityMap, AnnulusConvergesToTargetModuli)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    std::vector<Complex> annulus_centers = {Complex(0.0, 0.0)};
    std::vector<double> annulus_radii = {0.15};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        annulus_centers, annulus_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Annulus should converge within "
        << config.max_newton_iterations << " iterations";

    const auto& actual_radii = method.mp_canonical_domain->getHoleRadii();
    const auto& actual_centers = method.mp_canonical_domain->getHoleCenters();

    EXPECT_GT(actual_radii[0], 0.0) << "Converged inner radius must be positive";
    // For off-center annulus (inner at (0.3,0) r=0.15), the conformal radius is a
    // conformal invariant that depends on the eccentricity of the user domain.
    EXPECT_GT(actual_radii[0], 0.05) << "Converged radius should be physically reasonable";
    EXPECT_LT(actual_radii[0], 0.5) << "Converged radius should be physically reasonable";
    // Canonical annulus: inner circle centered at origin by convention (concentric representation)
    EXPECT_NEAR(std::abs(actual_centers[0]), 0.0, 1e-10)
        << "Annulus center should remain at origin";
}

// ============================================================================
// Boundary Correspondence Tests
// ============================================================================

class ConstructiveBoundaryCorrespondence : public ConstructiveValidationBase {};

TEST_F(ConstructiveBoundaryCorrespondence, IdentityM4MapPreservesBoundaries)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    std::vector<Complex> hole_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
    std::vector<double> hole_radii = {0.25, 0.25, 0.25};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        hole_centers, hole_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Identity m=4 must converge for boundary correspondence test";

    method.computeFourierCoefficients();

    // Evaluate map at points on the outer canonical boundary (unit circle)
    const int N_test = 32;
    const double outer_tol = 1e-4;
    for (int j = 0; j < N_test; ++j)
    {
        double theta = 2.0 * M_PI * j / N_test;
        Complex z = Complex(std::cos(theta), std::sin(theta));
        Complex w = method.map(z);
        double dist_from_unit_circle = std::abs(std::abs(w) - 1.0);
        EXPECT_LT(dist_from_unit_circle, outer_tol)
            << "Outer boundary point at theta=" << theta
            << ": |w|=" << std::abs(w) << ", expected ~1.0";
    }

    // Evaluate map at points on inner canonical boundaries
    const auto& converged_centers = method.mp_canonical_domain->getHoleCenters();
    const auto& converged_radii = method.mp_canonical_domain->getHoleRadii();

    const std::vector<Complex> target_centers = {Complex(-0.5, 0.0), Complex(0.25, 0.43), Complex(0.25, -0.43)};
    const double target_radius = 0.25;
    const double inner_tol = 1e-4;

    for (size_t nu = 0; nu < converged_centers.size(); ++nu)
    {
        for (int j = 0; j < N_test; ++j)
        {
            double theta = 2.0 * M_PI * j / N_test;
            Complex z = converged_centers[nu]
                        + converged_radii[nu] * Complex(std::cos(theta), std::sin(theta));
            Complex w = method.map(z);
            double dist_from_target = std::abs(std::abs(w - target_centers[nu]) - target_radius);
            EXPECT_LT(dist_from_target, inner_tol)
                << "Inner boundary " << nu << " at theta=" << theta
                << ": |w - c_target| = " << std::abs(w - target_centers[nu])
                << ", expected ~" << target_radius;
        }
    }
}

TEST_F(ConstructiveBoundaryCorrespondence, AnnulusMapPreservesBoundaries)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    std::vector<Complex> annulus_centers = {Complex(0.0, 0.0)};
    std::vector<double> annulus_radii = {0.15};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        annulus_centers, annulus_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Annulus must converge for boundary correspondence test";

    method.computeFourierCoefficients();

    // Points on outer canonical boundary (unit circle) should map to unit circle
    const int N_test = 32;
    const double tol = 1e-4;
    for (int j = 0; j < N_test; ++j)
    {
        double theta = 2.0 * M_PI * j / N_test;
        Complex z = Complex(std::cos(theta), std::sin(theta));
        Complex w = method.map(z);
        EXPECT_NEAR(std::abs(w), 1.0, tol)
            << "Outer boundary: |map(z)| at theta=" << theta;
    }

    // Points on inner canonical boundary should map to user inner boundary
    // (circle centered at 0.3+0i, radius 0.15)
    const auto& converged_radii = method.mp_canonical_domain->getHoleRadii();
    const auto& converged_centers = method.mp_canonical_domain->getHoleCenters();
    for (int j = 0; j < N_test; ++j)
    {
        double theta = 2.0 * M_PI * j / N_test;
        Complex z = converged_centers[0]
                    + converged_radii[0] * Complex(std::cos(theta), std::sin(theta));
        Complex w = method.map(z);
        double dist = std::abs(std::abs(w - Complex(0.3, 0)) - 0.15);
        EXPECT_LT(dist, tol)
            << "Inner boundary: |map(z) - center| at theta=" << theta
            << " is " << std::abs(w - Complex(0.3, 0)) << ", expected ~0.15";
    }
}

TEST_F(ConstructiveBoundaryCorrespondence, GeneralM3MapPreservesBoundaries)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );

    std::vector<Complex> m3_centers = {Complex(0.3, 0.2), Complex(-0.3, -0.1)};
    std::vector<double> m3_radii = {0.1, 0.12};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        m3_centers, m3_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "General m=3 must converge for boundary correspondence test";

    method.computeFourierCoefficients();

    // Points on outer canonical boundary (unit circle) should map to unit circle
    const int N_test = 32;
    const double tol = 1e-4;
    for (int j = 0; j < N_test; ++j)
    {
        double theta = 2.0 * M_PI * j / N_test;
        Complex z = Complex(std::cos(theta), std::sin(theta));
        Complex w = method.map(z);
        EXPECT_NEAR(std::abs(w), 1.0, tol)
            << "Outer boundary: |map(z)| at theta=" << theta;
    }

    // Points on inner canonical boundaries should map to corresponding user boundaries
    const auto& converged_centers = method.mp_canonical_domain->getHoleCenters();
    const auto& converged_radii = method.mp_canonical_domain->getHoleRadii();

    const std::vector<Complex> target_centers = {Complex(0.3, 0.2), Complex(-0.3, -0.1)};
    const std::vector<double> target_radii = {0.1, 0.12};

    for (size_t nu = 0; nu < converged_centers.size(); ++nu)
    {
        for (int j = 0; j < N_test; ++j)
        {
            double theta = 2.0 * M_PI * j / N_test;
            Complex z = converged_centers[nu]
                        + converged_radii[nu] * Complex(std::cos(theta), std::sin(theta));
            Complex w = method.map(z);
            double dist_from_target = std::abs(std::abs(w - target_centers[nu]) - target_radii[nu]);
            EXPECT_LT(dist_from_target, tol)
                << "Inner boundary " << nu << " at theta=" << theta
                << ": |w - c_target| = " << std::abs(w - target_centers[nu])
                << ", expected ~" << target_radii[nu];
        }
    }
}

// ============================================================================
// Conformal Moduli Consistency Tests
// ============================================================================

class ConstructiveModuliConsistency : public ConstructiveValidationBase {};

TEST_F(ConstructiveModuliConsistency, AnnulusModuliValid)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    std::vector<Complex> annulus_centers = {Complex(0.0, 0.0)};
    std::vector<double> annulus_radii = {0.15};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        annulus_centers, annulus_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged);

    const auto& radii = method.mp_canonical_domain->getHoleRadii();
    const auto& centers = method.mp_canonical_domain->getHoleCenters();

    for (size_t i = 0; i < radii.size(); ++i)
    {
        EXPECT_GT(radii[i], 0.0) << "Hole " << i << " radius must be positive";
    }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        double dist_plus_radius = std::abs(centers[i]) + radii[i];
        EXPECT_LT(dist_plus_radius, 1.0)
            << "Hole " << i << " must be inside unit disk: |c| + rho = " << dist_plus_radius;
    }
}

TEST_F(ConstructiveModuliConsistency, GeneralM3ModuliValid)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );

    std::vector<Complex> m3_centers = {Complex(0.3, 0.2), Complex(-0.3, -0.1)};
    std::vector<double> m3_radii = {0.1, 0.12};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        m3_centers, m3_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged);

    const auto& radii = method.mp_canonical_domain->getHoleRadii();
    const auto& centers = method.mp_canonical_domain->getHoleCenters();

    for (size_t i = 0; i < radii.size(); ++i)
    {
        EXPECT_GT(radii[i], 0.0) << "Hole " << i << " radius must be positive";
    }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        double dist_plus_radius = std::abs(centers[i]) + radii[i];
        EXPECT_LT(dist_plus_radius, 1.0)
            << "Hole " << i << " must be inside unit disk: |c| + rho = " << dist_plus_radius;
    }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        for (size_t j = i + 1; j < centers.size(); ++j)
        {
            double center_dist = std::abs(centers[i] - centers[j]);
            double radii_sum = radii[i] + radii[j];
            EXPECT_GT(center_dist, radii_sum)
                << "Holes " << i << " and " << j << " overlap: |c_i - c_j| = "
                << center_dist << " <= rho_i + rho_j = " << radii_sum;
        }
    }
}

TEST_F(ConstructiveModuliConsistency, IdentityM4ModuliValid)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    std::vector<Complex> hole_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
    std::vector<double> hole_radii = {0.25, 0.25, 0.25};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        hole_centers, hole_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged);

    const auto& radii = method.mp_canonical_domain->getHoleRadii();
    const auto& centers = method.mp_canonical_domain->getHoleCenters();

    for (size_t i = 0; i < radii.size(); ++i)
    {
        EXPECT_GT(radii[i], 0.0) << "Hole " << i << " radius must be positive";
    }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        double dist_plus_radius = std::abs(centers[i]) + radii[i];
        EXPECT_LT(dist_plus_radius, 1.0)
            << "Hole " << i << " must be inside unit disk: |c| + rho = " << dist_plus_radius;
    }

    for (size_t i = 0; i < centers.size(); ++i)
    {
        for (size_t j = i + 1; j < centers.size(); ++j)
        {
            double center_dist = std::abs(centers[i] - centers[j]);
            double radii_sum = radii[i] + radii[j];
            EXPECT_GT(center_dist, radii_sum)
                << "Holes " << i << " and " << j << " overlap: |c_i - c_j| = "
                << center_dist << " <= rho_i + rho_j = " << radii_sum;
        }
    }
}

// ============================================================================
// Convergence Rate Tests
// ============================================================================

class ConstructiveConvergenceRate : public ConstructiveValidationBase {};

TEST_F(ConstructiveConvergenceRate, AnnulusConvergesWithDecreasingResidual)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    std::vector<Complex> annulus_centers = {Complex(0.0, 0.0)};
    std::vector<double> annulus_radii = {0.15};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        annulus_centers, annulus_radii, config.N
    );
    method.initializeNewtonIteration();

    std::vector<double> residuals;
    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        residuals.push_back(method.getCurrentResidual());
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }

    ASSERT_TRUE(converged) << "Annulus should converge";
    ASSERT_GE(residuals.size(), 2u) << "Need at least 2 iterations to check trend";

    EXPECT_LT(residuals.back(), config.newton_tolerance);

    EXPECT_LT(residuals.back(), residuals.front())
        << "Final residual (" << residuals.back()
        << ") should be less than initial (" << residuals.front() << ")";
}

TEST_F(ConstructiveConvergenceRate, IdentityM4ConvergesWithDecreasingResidual)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    std::vector<Complex> hole_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
    std::vector<double> hole_radii = {0.25, 0.25, 0.25};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        hole_centers, hole_radii, config.N
    );
    method.initializeNewtonIteration();

    std::vector<double> residuals;
    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        residuals.push_back(method.getCurrentResidual());
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }

    ASSERT_TRUE(converged) << "Identity m=4 should converge";
    ASSERT_GE(residuals.size(), 2u) << "Need at least 2 iterations to check trend";

    EXPECT_LT(residuals.back(), config.newton_tolerance);

    EXPECT_LT(residuals.back(), residuals.front())
        << "Final residual (" << residuals.back()
        << ") should be less than initial (" << residuals.front() << ")";
}

TEST_F(ConstructiveConvergenceRate, GeneralM3ConvergesWithDecreasingResidual)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );

    std::vector<Complex> m3_centers = {Complex(0.3, 0.2), Complex(-0.3, -0.1)};
    std::vector<double> m3_radii = {0.1, 0.12};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        m3_centers, m3_radii, config.N
    );
    method.initializeNewtonIteration();

    std::vector<double> residuals;
    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        residuals.push_back(method.getCurrentResidual());
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }

    ASSERT_TRUE(converged) << "General m=3 should converge";

    EXPECT_LT(residuals.back(), config.newton_tolerance);

    // If multiple iterations occurred, verify overall decreasing trend
    if (residuals.size() >= 2)
    {
        EXPECT_LT(residuals.back(), residuals.front())
            << "Final residual (" << residuals.back()
            << ") should be less than initial (" << residuals.front() << ")";
    }
}

// ============================================================================
// Fourier Coefficient Properties Tests
// ============================================================================

class ConstructiveFourierProperties : public ConstructiveValidationBase {};

TEST_F(ConstructiveFourierProperties, IdentityM4OuterCoefficientIsUnity)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    std::vector<Complex> hole_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
    std::vector<double> hole_radii = {0.25, 0.25, 0.25};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        hole_centers, hole_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Identity m=4 must converge for Fourier test";

    method.computeFourierCoefficients();

    const auto& a = method.m_a;
    const int N = a.rows();

    // Outer boundary (column 0): the k=1 Fourier mode.
    // For the identity map on the unit circle, f(z) = z, so the
    // k=1 coefficient should have magnitude ~1.0 (the outer radius).
    double dominant_outer = std::abs(a(1, 0));
    EXPECT_NEAR(dominant_outer, 1.0, 0.01)
        << "Outer boundary k=1 Fourier coefficient: |a(1,0)| = " << dominant_outer;

    // Spectral decay: higher-order Fourier modes (indices N/4 to N/2-1)
    // should be smaller than the dominant k=1 mode
    double tail_max = 0.0;
    for (int j = N / 4; j < N / 2; ++j)
    {
        tail_max = std::max(tail_max, std::abs(a(j, 0)));
    }
    EXPECT_LT(tail_max, dominant_outer)
        << "High-order coefficients (max=" << tail_max
        << ") should be smaller than dominant (" << dominant_outer << ")";
}

TEST_F(ConstructiveFourierProperties, AnnulusFourierCoefficientDecay)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    std::vector<Complex> annulus_centers = {Complex(0.0, 0.0)};
    std::vector<double> annulus_radii = {0.15};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 2;
    method.m_is_annulus = true;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        annulus_centers, annulus_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "Annulus must converge for Fourier test";

    method.computeFourierCoefficients();

    const auto& a = method.m_a;
    const int N = a.rows();

    double dominant_outer = std::abs(a(1, 0));
    EXPECT_GT(dominant_outer, 0.1)
        << "Outer boundary should have significant dominant coefficient";

    // Spectral decay: high-frequency tail (indices 3N/8 to N/2-1) should be much smaller
    double tail_max = 0.0;
    for (int j = 3 * N / 8; j < N / 2; ++j)
    {
        tail_max = std::max(tail_max, std::abs(a(j, 0)));
    }
    EXPECT_LT(tail_max, 0.1 * dominant_outer)
        << "Tail coefficients (max=" << tail_max
        << ") should decay relative to dominant (" << dominant_outer << ")";
}

TEST_F(ConstructiveFourierProperties, GeneralM3FourierCoefficientDecay)
{
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );

    std::vector<Complex> m3_centers = {Complex(0.3, 0.2), Complex(-0.3, -0.1)};
    std::vector<double> m3_radii = {0.1, 0.12};

    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 3;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        m3_centers, m3_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();
        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    ASSERT_TRUE(converged) << "General m=3 must converge for Fourier test";

    method.computeFourierCoefficients();

    const auto& a = method.m_a;
    const int N = a.rows();
    const int m = a.cols();

    // Verify spectral decay for each boundary component
    for (int col = 0; col < m; ++col)
    {
        SCOPED_TRACE("Boundary component " + std::to_string(col));

        double dominant = std::abs(a(1, col));
        EXPECT_GT(dominant, 0.01)
            << "k=1 coefficient should be significant";

        // High-frequency tail should be much smaller than dominant mode
        double tail_max = 0.0;
        for (int j = 3 * N / 8; j < N / 2; ++j)
        {
            tail_max = std::max(tail_max, std::abs(a(j, col)));
        }
        EXPECT_LT(tail_max, 0.1 * dominant)
            << "Tail coefficients (max=" << tail_max
            << ") should decay relative to dominant (" << dominant << ")";
    }
}

// ============================================================================
// GH-92: Thesis3 N=256 Iteration Divergence Diagnostic
// ============================================================================

class Thesis3DiagnosticTest : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;
    std::shared_ptr<MultiplyConnectedDomain> domain;
    std::vector<Complex> initial_centers;
    std::vector<double> initial_radii;

    void SetUp() override
    {
        config.N = 256;
        config.max_newton_iterations = 20;
        config.newton_tolerance = 1e-14;
        config.cgm_tolerance = 1e-15;
        config.max_cgm_iterations = 100;
        config.verbose = false;

        auto outer = createCircularBoundary(Complex(0, 0), 1.0);
        auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
        auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
        auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
        domain = std::make_shared<MultiplyConnectedDomain>(
            std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
        );

        initial_centers = {Complex(-0.4, 0.0), Complex(0.35, 0.43), Complex(0.35, -0.43)};
        initial_radii = {0.25, 0.25, 0.25};
    }
};

TEST_F(Thesis3DiagnosticTest, Iteration1MatchesOctaveState)
{
    // Octave reference values (from thesis3_iteration_debug.m at N=256):
    // Iteration 1: normU = 0.1
    // After iteration 1: c = [-0.5; 0.25+0.43i; 0.25-0.43i], rho = [0.25; 0.25; 0.25]
    // (identity map for circles: first iteration corrects the offset initial guesses,
    //  second iteration confirms convergence at ~1e-14)
    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        initial_centers, initial_radii, config.N
    );
    method.initializeNewtonIteration();

    // --- Iteration 1 ---
    method.formSystem();

    // Verify D matrix and g vector are finite after iteration 1
    ASSERT_TRUE(method.m_D.allFinite())
        << "D matrix has non-finite values after iteration 1 formSystem";
    ASSERT_TRUE(method.m_g.allFinite())
        << "g vector has non-finite values after iteration 1 formSystem";

    method.solveSystem();

    // Verify solution is finite
    ASSERT_TRUE(method.m_U.allFinite())
        << "U vector has non-finite values after iteration 1 solveSystem";

    method.newtonUpdate();

    // Octave: normU for iteration 1 = 1.0e-01
    EXPECT_NEAR(method.m_current_residual, 0.1, 1e-3)
        << "Iteration 1 residual should be ~0.1 (matching Octave)";

    // After iteration 1, the moduli should be very close to target values
    // (identity map for circles: c -> actual centers, rho -> actual radii)
    const std::vector<Complex> target_centers = {
        Complex(-0.5, 0.0), Complex(0.25, 0.43), Complex(0.25, -0.43)};
    const std::vector<double> target_radii = {0.25, 0.25, 0.25};

    for (int i = 0; i < 3; ++i)
    {
        Complex c_val = method.m_conformal_moduli(2 * i);
        double rho_val = std::real(method.m_conformal_moduli(2 * i + 1));
        EXPECT_NEAR(std::abs(c_val - target_centers[i]), 0.0, 1e-6)
            << "After iter 1: c[" << i << "] = " << c_val
            << ", expected " << target_centers[i];
        EXPECT_NEAR(rho_val, target_radii[i], 1e-6)
            << "After iter 1: rho[" << i << "] = " << rho_val
            << ", expected " << target_radii[i];
    }

    // Verify S values are finite and in reasonable range
    ASSERT_TRUE(method.m_S.allFinite())
        << "S matrix has non-finite values after iteration 1";
}

TEST_F(Thesis3DiagnosticTest, Iteration2FormSystemProducesFiniteValues)
{
    // Regression test for GH-92: before the fix, iteration 2's D/g contained inf values
    // in C++. Octave reference: max|D2| = 256, max|g2| ~ 4.5e-14 (all finite).
    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        initial_centers, initial_radii, config.N
    );
    method.initializeNewtonIteration();

    // Run iteration 1
    method.formSystem();
    method.solveSystem();
    method.newtonUpdate();

    // --- Iteration 2 ---
    method.formSystem();

    // Verify D and g are finite (GH-92 regression point)
    double D_max = 0.0;
    for (int i = 0; i < method.m_D.rows(); ++i)
    {
        for (int j = 0; j < method.m_D.cols(); ++j)
        {
            D_max = std::max(D_max, std::abs(method.m_D(i, j)));
        }
    }
    double g_max = 0.0;
    for (int i = 0; i < method.m_g.size(); ++i)
    {
        g_max = std::max(g_max, std::abs(method.m_g(i)));
    }

    EXPECT_TRUE(method.m_D.allFinite())
        << "GH-92: D matrix has non-finite values in iteration 2. max|D| = " << D_max;
    EXPECT_TRUE(method.m_g.allFinite())
        << "GH-92: g vector has non-finite values in iteration 2. max|g| = " << g_max;

    // If D is finite, verify the solve also works
    if (method.m_D.allFinite() && method.m_g.allFinite())
    {
        method.solveSystem();
        EXPECT_TRUE(method.m_U.allFinite())
            << "GH-92: U vector has non-finite values in iteration 2 solve";

        method.newtonUpdate();

        // Octave: iteration 2 normU = 1.354e-14 (converged)
        EXPECT_LT(method.m_current_residual, 1e-10)
            << "Iteration 2 should have near-zero residual (Octave: 1.35e-14)";
    }
}

TEST_F(Thesis3DiagnosticTest, FullNewtonConvergenceN256)
{
    // Full convergence test at N=256 (thesis3 case that hung in CLI/GUI before GH-92 fix)
    FornbergMC method(config);
    method.setStatusManager(std::make_shared<StatusManager>());
    method.mp_user_domain = domain;
    method.m_connectivity = 4;
    method.m_is_annulus = false;
    method.mp_canonical_domain = std::make_shared<FornbergCanonicalDomain>(
        initial_centers, initial_radii, config.N
    );
    method.initializeNewtonIteration();

    bool converged = false;
    for (int iter = 0; iter < config.max_newton_iterations; ++iter)
    {
        method.formSystem();
        method.solveSystem();
        method.newtonUpdate();

        // Check for non-finite residual (the GH-92 failure mode)
        ASSERT_TRUE(std::isfinite(method.m_current_residual))
            << "GH-92: Non-finite residual at iteration " << (iter + 1)
            << ": " << method.m_current_residual;

        if (method.checkConvergence(config.newton_tolerance))
        {
            converged = true;
            break;
        }
    }
    EXPECT_TRUE(converged)
        << "Thesis3 (N=256) should converge within "
        << config.max_newton_iterations << " iterations";
}
