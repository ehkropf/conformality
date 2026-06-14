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
#include "../src/core/StatusManager.h"

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

// GH-52: Verify PMatrixBuilder rejects non-power-of-2 N values
TEST_F(PMatrixBuilderTest, ConstructorRejectsNonPowerOfTwoN)
{
    int connectivity = 3;
    bool is_annulus = false;
    FornbergMCConfiguration test_config = config;

    // Non-power-of-2 values should be rejected
    test_config.N = 7;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    test_config.N = 9;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    test_config.N = 15;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    test_config.N = 100;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    // Zero and negative values should also be rejected
    test_config.N = 0;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    test_config.N = -1;
    EXPECT_THROW(PMatrixBuilder(test_config, connectivity, is_annulus), std::invalid_argument);

    // Valid power-of-2 values should work
    // N=1 passes validation (2^0=1) but causes M=N/2=0, which is impractical
    test_config.N = 1;
    EXPECT_NO_THROW(PMatrixBuilder(test_config, connectivity, is_annulus));

    test_config.N = 2;  // Smallest practical value (M=1)
    EXPECT_NO_THROW(PMatrixBuilder(test_config, connectivity, is_annulus));

    test_config.N = 32;
    EXPECT_NO_THROW(PMatrixBuilder(test_config, connectivity, is_annulus));

    test_config.N = 128;
    EXPECT_NO_THROW(PMatrixBuilder(test_config, connectivity, is_annulus));
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
    // This was the bug from issue #34
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

// Regression test: Verify frequency indices are consistent across all components with
// connectivity > 2. The original bug fix was tested only with connectivity=2, which
// executes the initialization loop just twice. This test ensures the fix works across
// multiple loop iterations.
TEST_F(PMatrixBuilderTest, FrequencyIndicesMultipleComponents)
{
    config.N = 8;
    int connectivity = 4;  // Test with higher connectivity to verify loop executes multiple times
    PMatrixBuilder builder(config, connectivity, false);

    const auto& ref_indices = builder.getFrequencyIndices(0);

    // Verify all components have identical frequency ordering
    // This ensures the initialization loop in PMatrixBuilder::initializeFrequencyIndices()
    // works correctly for all connectivity values, not just connectivity=2
    for (int nu = 1; nu < connectivity; ++nu)
    {
        const auto& comp_indices = builder.getFrequencyIndices(nu);
        EXPECT_EQ(comp_indices, ref_indices)
            << "Component " << nu << " should have same frequency indices as component 0";

        // Also verify each component has correct FFTW ordering
        EXPECT_EQ(comp_indices[0], 0) << "Component " << nu << " should start with DC (0)";
        EXPECT_EQ(comp_indices[config.N / 2], -config.N / 2)
            << "Component " << nu << " Nyquist should be -N/2";
    }
}

// Regression test: Verify frequency index calculation handles N/2 integer division correctly
// for various power-of-2 values of N (FFTW requires N to be a power of 2). Tests that the
// Nyquist frequency bug fix works across different N values.
TEST_F(PMatrixBuilderTest, FrequencyIndicesPowerOfTwoN)
{
    int connectivity = 2;

    // Test with multiple power-of-2 values for N (up to maximum supported value)
    for (int N : {2, 4, 8, 16, 32, 64, 128, 256, 512})
    {
        config.N = N;
        PMatrixBuilder builder(config, connectivity, false);

        const auto& indices = builder.getFrequencyIndices(0);
        ASSERT_EQ(indices.size(), N) << "For N=" << N << ", size should be N";

        // Verify DC component
        EXPECT_EQ(indices[0], 0) << "For N=" << N << ", DC component should be 0";

        // Verify Nyquist frequency at index N/2 maps to -N/2
        EXPECT_EQ(indices[N / 2], -N / 2)
            << "For N=" << N << ", Nyquist at index " << N / 2 << " should be -" << N / 2;

        // Verify positive frequencies [0, N/2)
        for (int j = 0; j < N / 2; ++j)
        {
            EXPECT_EQ(indices[j], j)
                << "For N=" << N << ", positive frequency at index " << j << " should be " << j;
        }

        // Verify negative frequencies [N/2, N)
        for (int j = N / 2; j < N; ++j)
        {
            int expected = j - N;  // Maps to negative frequency
            EXPECT_EQ(indices[j], expected)
                << "For N=" << N << ", negative frequency at index " << j
                << " should be " << expected;
        }
    }
}

// Regression test: Verify annulus mode uses identical frequency indexing to general mode.
// While PMatrixBuilder::initializeFrequencyIndices() currently has no annulus-specific logic,
// this test prevents future regressions if annulus-specific optimizations are added.
TEST_F(PMatrixBuilderTest, FrequencyIndicesAnnulusConsistency)
{
    config.N = 8;
    PMatrixBuilder general_builder(config, 2, false);
    PMatrixBuilder annulus_builder(config, 2, true);

    const auto& general_indices = general_builder.getFrequencyIndices(0);
    const auto& annulus_indices = annulus_builder.getFrequencyIndices(0);

    EXPECT_EQ(annulus_indices, general_indices)
        << "Annulus mode should use identical frequency indexing to general mode";
    EXPECT_EQ(annulus_indices[config.N / 2], -config.N / 2)
        << "Annulus Nyquist frequency should be -N/2";
}

// Regression test for issue #35: Validates that the coincident circle center
// validation throws std::invalid_argument (changed from std::runtime_error to
// align with CLAUDE.md error handling guidelines for validation errors)
TEST_F(PMatrixBuilderTest, RejectsCoincidentCircleCentersGeneral)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Create moduli with two circles having nearly identical centers.
    // The separation of 1e-15 is below the 1e-14 tolerance threshold.
    ConformalModuli moduli;
    moduli.c.resize(2);
    moduli.rho.resize(2);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.c(1) = std::complex<double>(0.3, 0.2) + 1e-15;  // Below 1e-14 threshold
    moduli.rho(0) = 0.15;
    moduli.rho(1) = 0.12;

    EXPECT_THROW(builder.buildPMatrix(1, moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, RejectsCoincidentCircleCentersAnnulusNu1)
{
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, true);

    // Annulus mode baseline test: verifies annulus mode doesn't throw with valid input.
    // Note: With connectivity=2 and nu=1, the interaction loop in buildAnnulusPMatrix
    // has range [1, 0) and never executes, so the coincident center validation in that
    // interaction block is not tested by this case. This test serves as a regression
    // baseline to ensure annulus mode continues to work correctly.
    auto moduli = createAnnulusModuli(0.3);
    EXPECT_NO_THROW(builder.buildPMatrix(1, moduli));
}

TEST_F(PMatrixBuilderTest, RejectsCoincidentCircleCentersHigherConnectivity)
{
    int connectivity = 4;
    PMatrixBuilder builder(config, connectivity, false);

    // Test coincident center detection with higher connectivity (4 boundaries).
    // Circles at indices 1 and 2 are nearly coincident (below 1e-14 tolerance).
    // This exercises the validation in buildGeneralPMatrix for nu>=2 cases.
    ConformalModuli moduli;
    moduli.c.resize(3);
    moduli.rho.resize(3);
    moduli.c(0) = std::complex<double>(0.0, 0.0);
    moduli.c(1) = std::complex<double>(0.5, 0.3);
    moduli.c(2) = std::complex<double>(0.5, 0.3) + 1e-15;  // Below 1e-14 threshold
    moduli.rho(0) = 0.3;
    moduli.rho(1) = 0.15;
    moduli.rho(2) = 0.12;

    // Building P matrix for nu=2 triggers validation of c(1) vs c(2)
    EXPECT_THROW(builder.buildPMatrix(2, moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, RejectsCoincidentCircleCentersMultipleLoopIterations)
{
    int connectivity = 5;
    PMatrixBuilder builder(config, connectivity, false);

    // Test that validation works correctly across multiple loop iterations.
    // Circles at indices 1 and 3 are coincident, but 0 and 2 are not.
    // When building nu=2, loop tests L=0,1,3 (skips L=2=nu-1).
    // Should detect coincidence at L=3 iteration, not just first iteration.
    ConformalModuli moduli;
    moduli.c.resize(4);
    moduli.rho.resize(4);
    moduli.c(0) = std::complex<double>(0.2, 0.1);
    moduli.c(1) = std::complex<double>(0.5, 0.3);
    moduli.c(2) = std::complex<double>(-0.3, 0.4);
    moduli.c(3) = std::complex<double>(0.5, 0.3) + 1e-15;  // ≈ c(1), below 1e-14 threshold
    moduli.rho.setConstant(0.12);

    EXPECT_THROW(builder.buildPMatrix(2, moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, CoincidentCircleErrorMessageIncludesDiagnostics)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Verify that the enhanced error messages include diagnostic information:
    // actual distance value, threshold, and actionable guidance.
    ConformalModuli moduli;
    moduli.c.resize(2);
    moduli.rho.resize(2);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.c(1) = moduli.c(0) + 1e-15;
    moduli.rho(0) = 0.15;
    moduli.rho(1) = 0.12;

    try
    {
        builder.buildPMatrix(1, moduli);
        FAIL() << "Expected std::invalid_argument";
    }
    catch (const std::invalid_argument& e)
    {
        std::string msg = e.what();
        EXPECT_NE(msg.find("distance"), std::string::npos) << "Error message should include 'distance'";
        EXPECT_NE(msg.find("coincidence tolerance"), std::string::npos)
            << "Error message should reference the coincidence tolerance";
        EXPECT_NE(msg.find("Ensure circle centers are distinct"), std::string::npos)
            << "Error message should include actionable guidance";
    }
}

// Regression test for issue #36: Validates that ConformalModuli with incorrect size
// is rejected early with std::invalid_argument. For an m-connected domain, moduli
// must have exactly m-1 entries. Without validation, wrong-sized moduli cause
// out-of-bounds access (debug assertions or undefined behavior in release builds).
TEST_F(PMatrixBuilderTest, RejectsWrongSizedModuliTooSmall)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Create moduli with too few entries (1 instead of 2 for connectivity=3)
    ConformalModuli moduli;
    moduli.c.resize(1);
    moduli.rho.resize(1);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.rho(0) = 0.15;

    // Should throw std::invalid_argument for size mismatch
    EXPECT_THROW(builder.buildPMatrix(0, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildPMatrix(1, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildAllPMatrices(moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, RejectsWrongSizedModuliTooLarge)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Create moduli with too many entries (3 instead of 2 for connectivity=3)
    ConformalModuli moduli;
    moduli.c.resize(3);
    moduli.rho.resize(3);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.c(1) = std::complex<double>(-0.4, 0.3);
    moduli.c(2) = std::complex<double>(0.1, -0.3);
    moduli.rho(0) = 0.15;
    moduli.rho(1) = 0.12;
    moduli.rho(2) = 0.10;

    // Should throw std::invalid_argument for size mismatch
    EXPECT_THROW(builder.buildPMatrix(0, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildAllPMatrices(moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, RejectsMismatchedModuliArraySizes)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Create moduli where c and rho have different sizes
    ConformalModuli moduli;
    moduli.c.resize(2);
    moduli.rho.resize(1);  // Mismatch: c has 2, rho has 1
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.c(1) = std::complex<double>(-0.4, 0.3);
    moduli.rho(0) = 0.15;

    // Should throw std::invalid_argument for size mismatch
    EXPECT_THROW(builder.buildPMatrix(0, moduli), std::invalid_argument);
}

TEST_F(PMatrixBuilderTest, ModuliSizeValidationErrorMessage)
{
    int connectivity = 3;
    PMatrixBuilder builder(config, connectivity, false);

    // Verify error message includes expected size information
    ConformalModuli moduli;
    moduli.c.resize(1);
    moduli.rho.resize(1);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.rho(0) = 0.15;

    try
    {
        builder.buildPMatrix(0, moduli);
        FAIL() << "Expected std::invalid_argument";
    }
    catch (const std::invalid_argument& e)
    {
        std::string msg = e.what();
        EXPECT_NE(msg.find("moduli size mismatch"), std::string::npos)
            << "Error message should mention 'moduli size mismatch'";
        EXPECT_NE(msg.find("expected 2 entries"), std::string::npos)
            << "Error message should include expected size (connectivity-1)";
    }
}

// Critical test: Verify rejection of completely empty moduli arrays.
// For connectivity=2 (annulus), expected moduli size is 1, but empty arrays
// (size 0) should be rejected to prevent out-of-bounds access.
TEST_F(PMatrixBuilderTest, RejectsEmptyModuliForConnectivityTwo)
{
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, false);

    // Empty moduli arrays (should have size 1 for connectivity=2)
    ConformalModuli moduli;
    moduli.c.resize(0);
    moduli.rho.resize(0);

    EXPECT_THROW(builder.buildPMatrix(0, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildPMatrix(1, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildAllPMatrices(moduli), std::invalid_argument);
}

// Critical test: Verify all component indices reject wrong-sized moduli.
// Tests with higher connectivity (4 boundaries) to ensure the validation
// works correctly for all nu values, not just the first few.
TEST_F(PMatrixBuilderTest, RejectsWrongSizedModuliForAllComponents)
{
    int connectivity = 4;  // Higher connectivity for thorough testing
    PMatrixBuilder builder(config, connectivity, false);

    // Moduli too small (1 entry instead of 3)
    ConformalModuli moduli;
    moduli.c.resize(1);
    moduli.rho.resize(1);
    moduli.c(0) = std::complex<double>(0.3, 0.2);
    moduli.rho(0) = 0.15;

    // Verify each component index rejects wrong-sized moduli
    for (int nu = 0; nu < connectivity; ++nu)
    {
        EXPECT_THROW(builder.buildPMatrix(nu, moduli), std::invalid_argument)
            << "Component " << nu << " should reject wrong-sized moduli";
    }

    // Also verify buildAllPMatrices rejects
    EXPECT_THROW(builder.buildAllPMatrices(moduli), std::invalid_argument);
}

// Important test: Verify annulus mode rejects wrong-sized moduli.
// While validation happens before the mode split, this test ensures
// the validation path works correctly in annulus mode and prevents
// regressions from future refactoring.
TEST_F(PMatrixBuilderTest, RejectsWrongSizedModuliInAnnulusMode)
{
    int connectivity = 2;
    PMatrixBuilder builder(config, connectivity, true);  // Annulus mode

    // Empty moduli (should have size 1 for connectivity=2)
    ConformalModuli moduli;
    moduli.c.resize(0);
    moduli.rho.resize(0);

    EXPECT_THROW(builder.buildPMatrix(0, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildPMatrix(1, moduli), std::invalid_argument);
    EXPECT_THROW(builder.buildAllPMatrices(moduli), std::invalid_argument);
}

// =============================================================================
// StatusManager Integration Tests (GH-37)
// =============================================================================

class PMatrixBuilderStatusManagerTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        config.N = 64;
        config.max_newton_iterations = 10;
        config.max_cgm_iterations = 100;
        config.newton_tolerance = 1e-8;
        config.cgm_tolerance = 1e-8;
        config.verbose = false;
    }

    FornbergMCConfiguration config;
};

TEST_F(PMatrixBuilderStatusManagerTest, SetterGetterWorkCorrectly)
{
    auto statusManager = std::make_shared<StatusManager>();

    PMatrixBuilder builder(config, 3, false);

    // Initially has StrictNullStatusManager (never null)
    EXPECT_NE(builder.getStatusManager(), nullptr);
    EXPECT_EQ(std::dynamic_pointer_cast<StatusManager>(builder.getStatusManager()), nullptr)
        << "Default should not be a real StatusManager";

    // Set status manager
    builder.setStatusManager(statusManager);
    EXPECT_EQ(builder.getStatusManager(), statusManager);

    // Setting nullptr resets to StrictNullStatusManager (not null)
    builder.setStatusManager(nullptr);
    EXPECT_NE(builder.getStatusManager(), nullptr);
    EXPECT_EQ(std::dynamic_pointer_cast<StatusManager>(builder.getStatusManager()), nullptr)
        << "After clearing, should revert to StrictNullStatusManager";
}

TEST_F(PMatrixBuilderStatusManagerTest, LogsMessagesOnStateChange)
{
    auto statusManager = std::make_shared<StatusManager>();

    PMatrixBuilder builder(config, 3, false);
    builder.setStatusManager(statusManager);
    builder.setConnectivity(4);

    // Check that messages were logged
    const auto& messages = statusManager->getMessages();
    bool foundPMatrixBuilderMessage = false;
    for (const auto& msg : messages)
    {
        if (msg.component == "PMatrixBuilder")
        {
            foundPMatrixBuilderMessage = true;
            break;
        }
    }
    EXPECT_TRUE(foundPMatrixBuilderMessage) << "Expected messages with component 'PMatrixBuilder'";
}

TEST_F(PMatrixBuilderStatusManagerTest, NoExceptionWithoutStatusManager)
{
    // Verify all operations work without a StatusManager (null-safety)
    PMatrixBuilder builder(config, 3, false);

    // These should all work without exceptions
    EXPECT_NO_THROW(builder.setConnectivity(4));
    EXPECT_NO_THROW(builder.setAnnulusMode(false));

    // setAnnulusMode with rejected case (connectivity != 2) -- StrictNullStatusManager throws on WARNING
    EXPECT_THROW(builder.setAnnulusMode(true), std::runtime_error);
}
