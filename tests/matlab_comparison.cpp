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
#include "reference_data_loader.h"
#include "../src/methods/FornbergMC.h"
#include "../src/methods/FornbergMCConfiguration.h"
#include "../src/methods/PMatrixBuilder.h"
#include "../src/methods/ConformalModuli.h"
#include "../src/numerics/CGSolver.h"
#include "../src/domains/FornbergCanonicalDomain.h"
#include "../src/domains/Domain.h"

#include <cmath>
#include <string>

#ifndef REFERENCE_DATA_DIR
#error "REFERENCE_DATA_DIR must be defined"
#endif

static std::string refDataPath(const std::string& filename)
{
    return std::string(REFERENCE_DATA_DIR) + "/" + filename;
}

// Helper: compare two complex matrices element-wise with tolerance
static void expectComplexMatrixNear(const Eigen::MatrixXcd& actual,
                                    const Eigen::MatrixXcd& expected,
                                    double tol, const std::string& label)
{
    ASSERT_EQ(actual.rows(), expected.rows())
        << label << ": row count mismatch (" << actual.rows() << " vs " << expected.rows() << ")";
    ASSERT_EQ(actual.cols(), expected.cols())
        << label << ": col count mismatch (" << actual.cols() << " vs " << expected.cols() << ")";

    double max_diff = (actual - expected).cwiseAbs().maxCoeff();
    EXPECT_LT(max_diff, tol)
        << label << ": max element-wise difference = " << max_diff
        << " (tolerance = " << tol << ")"
        << "\n  actual norm = " << actual.norm()
        << "\n  expected norm = " << expected.norm();
}

// Helper: compare two complex vectors element-wise with tolerance
static void expectComplexVectorNear(const Eigen::VectorXcd& actual,
                                    const Eigen::VectorXcd& expected,
                                    double tol, const std::string& label)
{
    ASSERT_EQ(actual.size(), expected.size())
        << label << ": size mismatch (" << actual.size() << " vs " << expected.size() << ")";

    double max_diff = (actual - expected).cwiseAbs().maxCoeff();
    EXPECT_LT(max_diff, tol)
        << label << ": max element-wise difference = " << max_diff
        << " (tolerance = " << tol << ")"
        << "\n  actual norm = " << actual.norm()
        << "\n  expected norm = " << expected.norm()
        << "\n  norm ratio = " << (expected.norm() > 0 ? actual.norm() / expected.norm() : 0.0);
}

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
// Stage 1: P Matrix Comparison
// ============================================================================

class MatlabComparisonPMatrix : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;

    void SetUp() override
    {
        config.N = 64;
        config.verbose = false;
    }
};

TEST_F(MatlabComparisonPMatrix, AnnulusNu1)
{
    ReferenceDataLoader ref(refDataPath("annulus_P_nu1.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(1);
    moduli.c(0) = Complex(0.3, 0.0);
    moduli.rho = Eigen::VectorXd(1);
    moduli.rho(0) = 0.15;

    PMatrixBuilder builder(config, 2, true);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(0, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "Annulus P_nu1");
}

TEST_F(MatlabComparisonPMatrix, AnnulusNu2)
{
    ReferenceDataLoader ref(refDataPath("annulus_P_nu2.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(1);
    moduli.c(0) = Complex(0.3, 0.0);
    moduli.rho = Eigen::VectorXd(1);
    moduli.rho(0) = 0.15;

    PMatrixBuilder builder(config, 2, true);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(1, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "Annulus P_nu2");
}

TEST_F(MatlabComparisonPMatrix, GeneralM3Nu1)
{
    ReferenceDataLoader ref(refDataPath("general_m3_P_nu1.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(2);
    moduli.c(0) = Complex(0.3, 0.2);
    moduli.c(1) = Complex(-0.3, -0.1);
    moduli.rho = Eigen::VectorXd(2);
    moduli.rho(0) = 0.1;
    moduli.rho(1) = 0.12;

    PMatrixBuilder builder(config, 3, false);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(0, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "General m=3 P_nu1");
}

TEST_F(MatlabComparisonPMatrix, GeneralM3Nu2)
{
    ReferenceDataLoader ref(refDataPath("general_m3_P_nu2.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(2);
    moduli.c(0) = Complex(0.3, 0.2);
    moduli.c(1) = Complex(-0.3, -0.1);
    moduli.rho = Eigen::VectorXd(2);
    moduli.rho(0) = 0.1;
    moduli.rho(1) = 0.12;

    PMatrixBuilder builder(config, 3, false);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(1, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "General m=3 P_nu2");
}

TEST_F(MatlabComparisonPMatrix, GeneralM3Nu3)
{
    ReferenceDataLoader ref(refDataPath("general_m3_P_nu3.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(2);
    moduli.c(0) = Complex(0.3, 0.2);
    moduli.c(1) = Complex(-0.3, -0.1);
    moduli.rho = Eigen::VectorXd(2);
    moduli.rho(0) = 0.1;
    moduli.rho(1) = 0.12;

    PMatrixBuilder builder(config, 3, false);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(2, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "General m=3 P_nu3");
}

TEST_F(MatlabComparisonPMatrix, IdentityM4Nu1)
{
    ReferenceDataLoader ref(refDataPath("identity_m4_P_nu1.json"));
    Eigen::MatrixXcd P_expected = ref.getComplexMatrix("P");

    ConformalModuli moduli;
    moduli.c = Eigen::VectorXcd(3);
    moduli.c(0) = Complex(-0.4, 0.0);
    moduli.c(1) = Complex(0.35, 0.43);
    moduli.c(2) = Complex(0.35, -0.43);
    moduli.rho = Eigen::VectorXd(3);
    moduli.rho(0) = 0.25;
    moduli.rho(1) = 0.25;
    moduli.rho(2) = 0.25;

    PMatrixBuilder builder(config, 4, false);
    Eigen::MatrixXcd P_actual = builder.buildPMatrix(0, moduli);

    expectComplexMatrixNear(P_actual, P_expected, 1e-12, "Identity m=4 P_nu1");
}

// ============================================================================
// Stage 2: formSystem Comparison
// ============================================================================

class MatlabComparisonFormSystem : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;

    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-14;
        config.cgm_tolerance = 1e-15;
        config.max_cgm_iterations = 20;
        config.verbose = false;
    }
};

TEST_F(MatlabComparisonFormSystem, AnnulusDAndG)
{
    ReferenceDataLoader ref(refDataPath("annulus_form_system_iter1.json"));
    Eigen::MatrixXcd D_expected = ref.getComplexMatrix("D");
    Eigen::VectorXcd g_expected = ref.getComplexVector("g");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    // MATLAB annulus convention: c(1) = 0 (inner canonical circle centered at origin)
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
    method.formSystem();

    expectComplexMatrixNear(method.m_D, D_expected, 1e-10, "Annulus D matrix");
    expectComplexVectorNear(method.m_g, g_expected, 1e-10, "Annulus g vector");
}

TEST_F(MatlabComparisonFormSystem, GeneralM3DAndG)
{
    ReferenceDataLoader ref(refDataPath("general_m3_form_system_iter1.json"));
    Eigen::MatrixXcd D_expected = ref.getComplexMatrix("D");
    Eigen::VectorXcd g_expected = ref.getComplexVector("g");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(0.3, 0.2), 0.1);
    auto inner2 = createCircularBoundary(Complex(-0.3, -0.1), 0.12);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2}
    );

    // Use explicit canonical domain matching MATLAB initial guesses
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
    method.formSystem();

    expectComplexMatrixNear(method.m_D, D_expected, 1e-10, "General m=3 D matrix");
    expectComplexVectorNear(method.m_g, g_expected, 1e-10, "General m=3 g vector");
}

TEST_F(MatlabComparisonFormSystem, IdentityM4DAndG)
{
    ReferenceDataLoader ref(refDataPath("identity_m4_form_system_iter1.json"));
    Eigen::MatrixXcd D_expected = ref.getComplexMatrix("D");
    Eigen::VectorXcd g_expected = ref.getComplexVector("g");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner1 = createCircularBoundary(Complex(-0.5, 0.0), 0.25);
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3}
    );

    // Use th_gen_ex3 initial guesses: [-.4, .25], [.35+.43i, .25], [.35-.43i, .25]
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
    method.formSystem();

    expectComplexMatrixNear(method.m_D, D_expected, 1e-10, "Identity m=4 D matrix");
    expectComplexVectorNear(method.m_g, g_expected, 1e-10, "Identity m=4 g vector");
}

// ============================================================================
// Stage 3: solveSystem Comparison
// ============================================================================

class MatlabComparisonSolveSystem : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;

    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-14;
        config.cgm_tolerance = 1e-15;
        config.max_cgm_iterations = 20;
        config.verbose = false;
    }
};

TEST_F(MatlabComparisonSolveSystem, AnnulusU)
{
    ReferenceDataLoader ref(refDataPath("annulus_solve_system_iter1.json"));
    Eigen::VectorXcd U_expected = ref.getComplexVector("U");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    // MATLAB annulus convention: c(1) = 0 (inner canonical circle centered at origin)
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
    method.formSystem();
    method.solveSystem();

    Eigen::VectorXcd U_actual = method.m_U;

    // Report norm ratio to detect FFT normalization mismatches (manifests as factor-of-2)
    double norm_ratio = U_actual.norm() / U_expected.norm();
    EXPECT_NEAR(norm_ratio, 1.0, 0.1)
        << "Norm ratio (C++/MATLAB) = " << norm_ratio
        << " (if ~2.0, factor-of-2 discrepancy)";

    expectComplexVectorNear(U_actual, U_expected, 1e-6, "Annulus U vector");
}

TEST_F(MatlabComparisonSolveSystem, IdentityM4U)
{
    ReferenceDataLoader ref(refDataPath("identity_m4_solve_system_iter1.json"));
    Eigen::VectorXcd U_expected = ref.getComplexVector("U");

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
    method.formSystem();
    method.solveSystem();

    Eigen::VectorXcd U_actual = method.m_U;

    double norm_ratio = (U_expected.norm() > 0) ?
        U_actual.norm() / U_expected.norm() : U_actual.norm();
    if (U_expected.norm() > 1e-14)
    {
        EXPECT_NEAR(norm_ratio, 1.0, 0.1)
            << "Norm ratio (C++/MATLAB) = " << norm_ratio;
    }

    expectComplexVectorNear(U_actual, U_expected, 1e-6, "Identity m=4 U vector");
}

// ============================================================================
// Stage 4: newtonUpdate Comparison
// ============================================================================

class MatlabComparisonNewtonUpdate : public ::testing::Test
{
protected:
    FornbergMCConfiguration config;

    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.newton_tolerance = 1e-14;
        config.cgm_tolerance = 1e-15;
        config.max_cgm_iterations = 20;
        config.verbose = false;
    }
};

TEST_F(MatlabComparisonNewtonUpdate, AnnulusFirstIteration)
{
    ReferenceDataLoader ref(refDataPath("annulus_newton_update_iter1.json"));
    Eigen::MatrixXcd S_expected_cx = ref.getComplexMatrix("S");
    Eigen::MatrixXd S_expected = S_expected_cx.real();
    Eigen::VectorXcd c_expected = ref.getComplexVector("c");
    Eigen::VectorXd rho_expected = ref.getRealVectorFromComplex("rho");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    // MATLAB annulus convention: c(1) = 0 (inner canonical circle centered at origin)
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
    method.formSystem();
    method.solveSystem();
    method.newtonUpdate();

    // Compare updated S
    double S_max_diff = (method.m_S - S_expected).cwiseAbs().maxCoeff();
    EXPECT_LT(S_max_diff, 1e-6)
        << "Annulus S after newton update: max diff = " << S_max_diff;

    // Compare updated rho (from canonical domain)
    const auto& hole_radii = method.mp_canonical_domain->getHoleRadii();
    for (int i = 0; i < rho_expected.size(); ++i)
    {
        EXPECT_NEAR(hole_radii[i], rho_expected(i), 1e-6)
            << "rho[" << i << "] mismatch";
    }

    // Compare updated c (from canonical domain)
    const auto& hole_centers = method.mp_canonical_domain->getHoleCenters();
    for (int i = 0; i < c_expected.size(); ++i)
    {
        EXPECT_NEAR(std::abs(hole_centers[i] - c_expected(i)), 0.0, 1e-6)
            << "c[" << i << "] mismatch: got " << hole_centers[i]
            << " expected " << c_expected(i);
    }
}

TEST_F(MatlabComparisonNewtonUpdate, IdentityM4FirstIteration)
{
    ReferenceDataLoader ref(refDataPath("identity_m4_newton_update_iter1.json"));
    Eigen::MatrixXcd S_expected_cx = ref.getComplexMatrix("S");
    Eigen::MatrixXd S_expected = S_expected_cx.real();
    Eigen::VectorXcd c_expected = ref.getComplexVector("c");
    Eigen::VectorXd rho_expected = ref.getRealVectorFromComplex("rho");

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
    method.formSystem();
    method.solveSystem();
    method.newtonUpdate();

    double S_max_diff = (method.m_S - S_expected).cwiseAbs().maxCoeff();
    EXPECT_LT(S_max_diff, 1e-6)
        << "Identity m=4 S after newton update: max diff = " << S_max_diff;

    const auto& actual_radii = method.mp_canonical_domain->getHoleRadii();
    for (int i = 0; i < rho_expected.size(); ++i)
    {
        EXPECT_NEAR(actual_radii[i], rho_expected(i), 1e-6)
            << "rho[" << i << "] mismatch";
    }

    const auto& actual_centers = method.mp_canonical_domain->getHoleCenters();
    for (int i = 0; i < c_expected.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual_centers[i] - c_expected(i)), 0.0, 1e-6)
            << "c[" << i << "] mismatch: got " << actual_centers[i]
            << " expected " << c_expected(i);
    }
}

// ============================================================================
// Stage 5: Full Convergence Comparison
// ============================================================================

class MatlabComparisonConvergence : public ::testing::Test
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

TEST_F(MatlabComparisonConvergence, AnnulusConverges)
{
    ReferenceDataLoader ref(refDataPath("annulus_converged.json"));
    Eigen::VectorXcd c_expected = ref.getComplexVector("c");
    Eigen::VectorXd rho_expected = ref.getRealVectorFromComplex("rho");
    Eigen::MatrixXcd a_expected = ref.getComplexMatrix("a");

    auto outer = createCircularBoundary(Complex(0, 0), 1.0);
    auto inner = createCircularBoundary(Complex(0.3, 0), 0.15);
    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner}
    );

    // MATLAB annulus convention: c(1) = 0 (inner canonical circle centered at origin)
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

    // Run full Newton iteration
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

    method.computeFourierCoefficients();

    // Compare conformal moduli
    const auto& actual_radii = method.mp_canonical_domain->getHoleRadii();
    for (int i = 0; i < rho_expected.size(); ++i)
    {
        EXPECT_NEAR(actual_radii[i], rho_expected(i), 1e-8)
            << "Final rho[" << i << "] mismatch";
    }
    const auto& actual_centers = method.mp_canonical_domain->getHoleCenters();
    for (int i = 0; i < c_expected.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual_centers[i] - c_expected(i)), 0.0, 1e-8)
            << "Final c[" << i << "] mismatch";
    }

    // Compare Fourier coefficients
    expectComplexMatrixNear(method.m_a, a_expected, 1e-8, "Annulus Fourier coefficients");
}

TEST_F(MatlabComparisonConvergence, IdentityM4Converges)
{
    ReferenceDataLoader ref(refDataPath("identity_m4_converged.json"));
    Eigen::VectorXcd c_expected = ref.getComplexVector("c");
    Eigen::VectorXd rho_expected = ref.getRealVectorFromComplex("rho");
    Eigen::MatrixXcd a_expected = ref.getComplexMatrix("a");

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

    method.computeFourierCoefficients();

    // For identity map (domain = canonical domain), converged c and rho should
    // match target domain boundary centers and radii within iteration tolerance
    const auto& actual_radii = method.mp_canonical_domain->getHoleRadii();
    for (int i = 0; i < rho_expected.size(); ++i)
    {
        EXPECT_NEAR(actual_radii[i], rho_expected(i), 1e-4)
            << "Final rho[" << i << "] mismatch";
    }
    const auto& actual_centers = method.mp_canonical_domain->getHoleCenters();
    for (int i = 0; i < c_expected.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual_centers[i] - c_expected(i)), 0.0, 1e-4)
            << "Final c[" << i << "] mismatch";
    }

    // Compare Fourier coefficients
    expectComplexMatrixNear(method.m_a, a_expected, 1e-4, "Identity m=4 Fourier coefficients");
}
