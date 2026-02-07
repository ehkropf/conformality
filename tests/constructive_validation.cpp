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

#include <cmath>

// Helper: create circular boundary
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
// Test 1: Identity Map Tests
// ============================================================================

class ConstructiveIdentityMap : public ::testing::Test
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

TEST_F(ConstructiveIdentityMap, IdentityM4ConvergesToTargetModuli)
{
    // th_gen_ex3 domain: 4-connected, circular boundaries
    // When target domain IS circular, the map is the identity and
    // converged moduli must match target boundary geometry
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
    EXPECT_NEAR(std::abs(actual_centers[0]), 0.0, 1e-10)
        << "Annulus center should remain at origin";
}

// ============================================================================
// Test 2: Boundary Correspondence Tests
// ============================================================================

class ConstructiveBoundaryCorrespondence : public ::testing::Test
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
    const double outer_tol = 5e-3;
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
    const double inner_tol = 5e-3;

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
    const double tol = 5e-3;
    for (int j = 0; j < N_test; ++j)
    {
        double theta = 2.0 * M_PI * j / N_test;
        Complex z = Complex(std::cos(theta), std::sin(theta));
        Complex w = method.map(z);
        EXPECT_NEAR(std::abs(w), 1.0, tol)
            << "Outer boundary: |map(z)| at theta=" << theta;
    }

    // Points on inner canonical boundary should map to the user inner boundary
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

// ============================================================================
// Test 3: Conformal Moduli Consistency Tests
// ============================================================================

class ConstructiveModuliConsistency : public ::testing::Test
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
// Test 4: Convergence Rate Tests
// ============================================================================

class ConstructiveConvergenceRate : public ::testing::Test
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
// Test 5: Fourier Coefficient Properties Tests
// ============================================================================

class ConstructiveFourierProperties : public ::testing::Test
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

TEST_F(ConstructiveFourierProperties, IdentityM4DominantCoefficientMatchesRadius)
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

    // Outer boundary (column 0): a(1,0) is the z^1 coefficient
    // For identity map of unit circle, this should have magnitude ~1.0
    double dominant_outer = std::abs(a(1, 0));
    EXPECT_NEAR(dominant_outer, 1.0, 0.1)
        << "Outer boundary dominant Fourier coefficient: |a(1,0)| = " << dominant_outer;

    // Spectral decay: higher-order coefficients should be smaller than dominant
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

    // Spectral decay: last quarter of coefficients should be much smaller
    double tail_max = 0.0;
    for (int j = 3 * N / 8; j < N / 2; ++j)
    {
        tail_max = std::max(tail_max, std::abs(a(j, 0)));
    }
    EXPECT_LT(tail_max, 0.1 * dominant_outer)
        << "Tail coefficients (max=" << tail_max
        << ") should decay relative to dominant (" << dominant_outer << ")";
}
