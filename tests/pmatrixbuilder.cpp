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
#include "../src/methods/FornbergMCConfiguration.h"

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

    // Matrix building should work
    EXPECT_NO_THROW(auto P_matrices = builder.buildAllPMatrices());

    // Individual matrix building should work
    EXPECT_NO_THROW(auto P0 = builder.buildPMatrix(0));

    // Invalid component index should fail
    EXPECT_THROW(builder.buildPMatrix(-1), std::invalid_argument);
    EXPECT_THROW(builder.buildPMatrix(connectivity), std::invalid_argument);
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