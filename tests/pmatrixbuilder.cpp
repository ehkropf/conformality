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
#include "../src/methods/PMatrixBuilder.h"
#include "../src/methods/ConformalModuli.h"
#include "../src/methods/FornbergMCConfiguration.h"

namespace
{
    constexpr double kMatrixElementTol = 1e-12;
    constexpr double kBlockComparisonTol = 1e-12;
}

class PMatrixBuilderTest : public ::testing::Test
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

    /// Create dummy conformal moduli for testing
    ConformalModuli createTestModuli(int connectivity)
    {
        ConformalModuli moduli;
        moduli.c.resize(connectivity - 1);
        moduli.rho.resize(connectivity - 1);

        // Simple test values: inner circles at c=0.3, -0.4+0.2i, etc. with radii 0.15, 0.12, etc.
        for (int i = 0; i < connectivity - 1; ++i)
        {
            moduli.c(i) = std::complex<double>(0.3 - 0.7 * i, 0.2 * i);
            moduli.rho(i) = 0.15 - 0.03 * i;
        }
        return moduli;
    }

    /// Create annulus moduli (single inner circle at origin)
    ConformalModuli createAnnulusModuli(double rho = 0.3)
    {
        ConformalModuli moduli;
        moduli.c.resize(1);
        moduli.rho.resize(1);
        moduli.c(0) = 0.0;  // Centered at origin for true annulus
        moduli.rho(0) = rho;
        return moduli;
    }

    FornbergMCConfiguration config;
};

// Test PMatrixBuilder
TEST_F(PMatrixBuilderTest, Construction)
{
    int connectivity = 3;
    bool is_annulus = false;

    // Construction should work
    EXPECT_NO_THROW(PMatrixBuilder builder(config, connectivity, is_annulus));

    PMatrixBuilder builder(config, connectivity, is_annulus);

    // Basic property checks
    EXPECT_EQ(builder.isAnnulusMode(), is_annulus);
    EXPECT_GT(builder.getSystemSize(), 0);
}

TEST_F(PMatrixBuilderTest, MatrixBuilding)
{
    int connectivity = 3;
    bool is_annulus = false;
    PMatrixBuilder builder(config, connectivity, is_annulus);
    auto moduli = createTestModuli(connectivity);

    // Matrix building should work
    EXPECT_NO_THROW(auto P_matrices = builder.buildAllPMatrices(moduli));

    // Individual matrix building should work
    EXPECT_NO_THROW(auto P0 = builder.buildPMatrix(0, moduli));

    // Invalid component index should fail
    EXPECT_THROW(builder.buildPMatrix(-1, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildPMatrix(connectivity, moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, AnnulusMode)
{
    int connectivity = 2; // Annulus case
    bool is_annulus = true;

    PMatrixBuilder builder(config, connectivity, is_annulus);

    EXPECT_TRUE(builder.isAnnulusMode());
    
    // Annulus should have different system size
    int annulus_size = builder.getSystemSize();
    
    // Compare with general case
    PMatrixBuilder general_builder(config, connectivity, false);
    int general_size = general_builder.getSystemSize();
    
    // Sizes may differ due to optimization
    EXPECT_GT(annulus_size, 0);
    EXPECT_GT(general_size, 0);
}

TEST_F(PMatrixBuilderTest, ConnectivityUpdate)
{
    int initial_connectivity = 2;
    PMatrixBuilder builder(config, initial_connectivity, false);
    
    int initial_size = builder.getSystemSize();
    
    // Update connectivity
    int new_connectivity = 3;
    EXPECT_NO_THROW(builder.setConnectivity(new_connectivity));
    
    int new_size = builder.getSystemSize();
    
    // System size should change with connectivity
    EXPECT_NE(initial_size, new_size);
}

TEST_F(PMatrixBuilderTest, AnnulusModeToggle)
{
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, false);

    EXPECT_FALSE(builder.isAnnulusMode());

    // Switch to annulus mode
    EXPECT_NO_THROW(builder.setAnnulusMode(true));
    EXPECT_TRUE(builder.isAnnulusMode());

    // Switch back
    EXPECT_NO_THROW(builder.setAnnulusMode(false));
    EXPECT_FALSE(builder.isAnnulusMode());
}

TEST_F(PMatrixBuilderTest, GeneralMatrixDimensions)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);
    auto moduli = createTestModuli(connectivity);

    int M = config.N / 2;
    int expected_rows = connectivity * M;
    int expected_cols = config.N;

    auto P_matrices = builder.buildAllPMatrices(moduli);

    EXPECT_EQ(P_matrices.size(), static_cast<size_t>(connectivity));
    for (int nu = 0; nu < connectivity; ++nu)
    {
        EXPECT_EQ(P_matrices[nu].rows(), expected_rows);
        EXPECT_EQ(P_matrices[nu].cols(), expected_cols);
    }
}

TEST_F(PMatrixBuilderTest, AnnulusMatrixDimensions)
{
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, true);
    auto moduli = createAnnulusModuli();

    int M = config.N / 2;
    int expected_rows = connectivity * M;
    int expected_cols = config.N;

    auto P_matrices = builder.buildAllPMatrices(moduli);

    EXPECT_EQ(P_matrices.size(), static_cast<size_t>(connectivity));
    for (int nu = 0; nu < connectivity; ++nu)
    {
        EXPECT_EQ(P_matrices[nu].rows(), expected_rows);
        EXPECT_EQ(P_matrices[nu].cols(), expected_cols);
    }
}

TEST_F(PMatrixBuilderTest, GeneralP0IdentityBlock)
{
    // For nu=0 (outer boundary), the first M rows, columns M to N-1
    // should contain an identity matrix
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);
    auto moduli = createTestModuli(connectivity);

    auto P0 = builder.buildPMatrix(0, moduli);

    int M = config.N / 2;
    Eigen::MatrixXcd expected_identity = Eigen::MatrixXcd::Identity(M, M);
    Eigen::MatrixXcd identity_block = P0.block(0, M, M, M);

    EXPECT_TRUE(identity_block.isApprox(expected_identity, kBlockComparisonTol))
        << "P0 identity block mismatch";
}

TEST_F(PMatrixBuilderTest, GeneralP1NegIdentityBlock)
{
    // For nu>=1 (inner boundaries), there should be a -Identity block
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);
    auto moduli = createTestModuli(connectivity);

    auto P1 = builder.buildPMatrix(1, moduli);

    int M = config.N / 2;
    Eigen::MatrixXcd expected_neg_identity = -Eigen::MatrixXcd::Identity(M, M);
    // For nu=1, the -Identity block is at rows M to 2M-1, columns 0 to M-1
    Eigen::MatrixXcd neg_identity_block = P1.block(M, 0, M, M);

    EXPECT_TRUE(neg_identity_block.isApprox(expected_neg_identity, kBlockComparisonTol))
        << "P1 negative identity block mismatch";
}

TEST_F(PMatrixBuilderTest, AnnulusWorkedExample)
{
    // Test case from plan22.md Appendix E: concentric annulus
    // N=8, M=4, rho=0.3
    config.N = 8;
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, true);
    auto moduli = createAnnulusModuli(0.3);

    auto P0 = builder.buildPMatrix(0, moduli);
    auto P1 = builder.buildPMatrix(1, moduli);

    int M = 4;
    double rho = 0.3;

    // Verify P0 dimensions
    EXPECT_EQ(P0.rows(), 8);
    EXPECT_EQ(P0.cols(), 8);

    // Verify P0 identity block in upper-right
    Eigen::MatrixXcd expected_identity = Eigen::MatrixXcd::Identity(M, M);
    EXPECT_TRUE(P0.block(0, M, M, M).isApprox(expected_identity, kBlockComparisonTol));

    // Verify P0 diagonal block (rows M to 2M-1, cols 0 to M-1)
    // Should be diag(rho^0, rho^1, rho^2, rho^3)
    for (int k = 0; k < M; ++k)
    {
        EXPECT_NEAR(P0(M + k, k).real(), std::pow(rho, k), kMatrixElementTol);
        EXPECT_NEAR(P0(M + k, k).imag(), 0.0, kMatrixElementTol);
    }

    // Verify P1 dimensions
    EXPECT_EQ(P1.rows(), 8);
    EXPECT_EQ(P1.cols(), 8);

    // Verify P1 -Identity block in lower-left
    EXPECT_TRUE(P1.block(M, 0, M, M).isApprox(-expected_identity, kBlockComparisonTol));

    // Verify P1 diagonal block (rows 0 to M-1, cols M to N-1)
    // Should be -diag(rho^M, rho^(M-1), ..., rho^1) = -diag(rho^4, rho^3, rho^2, rho^1)
    for (int k = 0; k < M; ++k)
    {
        double expected_val = -std::pow(rho, M - k);  // M=4, so rho^4, rho^3, rho^2, rho^1
        EXPECT_NEAR(P1(k, M + k).real(), expected_val, kMatrixElementTol);
        EXPECT_NEAR(P1(k, M + k).imag(), 0.0, kMatrixElementTol);
    }
}

TEST_F(PMatrixBuilderTest, NormalizationConditionsGeneral)
{
    // Test applyNormalizationConditions for general (non-annulus) case
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    int M = config.N / 2;
    int rows = connectivity * M + 2;  // Extra rows for normalization
    int cols = connectivity * config.N + 3 * (connectivity - 1);

    Eigen::MatrixXcd D = Eigen::MatrixXcd::Random(rows, cols);
    Eigen::VectorXcd g = Eigen::VectorXcd::Random(rows);
    double norm_cond_value = 1.5;

    builder.applyNormalizationConditions(D, g, norm_cond_value);

    // Last row should be [1, 0, 0, ..., 0]
    int last_row = rows - 1;
    EXPECT_NEAR(D(last_row, 0).real(), 1.0, kMatrixElementTol);
    for (int j = 1; j < cols; ++j)
    {
        EXPECT_NEAR(std::abs(D(last_row, j)), 0.0, kMatrixElementTol);
    }

    // Last element of g should be 0
    EXPECT_NEAR(std::abs(g(last_row)), 0.0, kMatrixElementTol);
}

TEST_F(PMatrixBuilderTest, NormalizationConditionsAnnulus)
{
    // For annulus case, applyNormalizationConditions should do nothing
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, true);

    int M = config.N / 2;
    int rows = connectivity * M;
    int cols = connectivity * config.N;

    Eigen::MatrixXcd D = Eigen::MatrixXcd::Random(rows, cols);
    Eigen::VectorXcd g = Eigen::VectorXcd::Random(rows);
    Eigen::MatrixXcd D_orig = D;
    Eigen::VectorXcd g_orig = g;

    builder.applyNormalizationConditions(D, g, 1.5);

    // Matrices should be unchanged
    EXPECT_TRUE(D.isApprox(D_orig, kMatrixElementTol));
    EXPECT_TRUE(g.isApprox(g_orig, kMatrixElementTol));
}

TEST_F(PMatrixBuilderTest, FrequencyIndicesFFTWOrder)
{
    // Test that frequency indices follow FFTW output order:
    // [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
    // This is a regression test for issue #34

    config.N = 8;  // Use N=8 for clear test case (N/2 = 4)
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, false);

    // Get frequency indices for first component
    const auto& freq_indices = builder.getFrequencyIndices(0);

    ASSERT_EQ(freq_indices.size(), config.N);

    // Positive frequencies and zero: indices 0 to N/2-1 should map to [0, 1, 2, 3]
    for (int j = 0; j < config.N / 2; ++j)
    {
        EXPECT_EQ(freq_indices[j], j)
            << "Positive frequency at index " << j << " should equal " << j;
    }

    // Negative frequencies including Nyquist: indices N/2 to N-1 should map to [-N/2, -N/2+1, ..., -1]
    // For N=8: indices 4,5,6,7 should map to [-4, -3, -2, -1]
    for (int j = config.N / 2; j < config.N; ++j)
    {
        int expected_freq = j - config.N;
        EXPECT_EQ(freq_indices[j], expected_freq)
            << "Negative frequency at index " << j << " should equal " << expected_freq;
    }

    // Critical test: Nyquist frequency at index N/2 should be -N/2, not N/2
    // This is the bug from issue #34
    EXPECT_EQ(freq_indices[config.N / 2], -config.N / 2)
        << "Nyquist frequency at index N/2=" << config.N / 2
        << " should be -N/2=" << -config.N / 2 << ", not N/2";

    // Verify all connectivity components have consistent frequency ordering
    for (int nu = 1; nu < connectivity; ++nu)
    {
        const auto& comp_indices = builder.getFrequencyIndices(nu);
        EXPECT_EQ(comp_indices, freq_indices)
            << "Component " << nu << " should have same frequency indices as component 0";
    }
}

