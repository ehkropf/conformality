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
#include "../src/examples/ThesisExamples.h"

using namespace conformality::examples;

// --- availableExamples ---

TEST(ThesisExamples, AvailableExamplesReturnsFourItems)
{
    auto examples = ThesisExamples::availableExamples();
    EXPECT_EQ(examples.size(), 4u);
}

TEST(ThesisExamples, AvailableExamplesOrderedSimplestFirst)
{
    auto examples = ThesisExamples::availableExamples();
    ASSERT_EQ(examples.size(), 4u);
    EXPECT_EQ(examples[0], 3);  // Identity (simplest)
    EXPECT_EQ(examples[1], 5);  // Ellipses
    EXPECT_EQ(examples[2], 2);  // Mixed (inverted ellipse)
    EXPECT_EQ(examples[3], 4);  // High connectivity
}

// --- Invalid example number ---

TEST(ThesisExamples, InvalidExampleThrows)
{
    EXPECT_THROW(ThesisExamples::getExample(1), std::invalid_argument);
    EXPECT_THROW(ThesisExamples::getExample(0), std::invalid_argument);
    EXPECT_THROW(ThesisExamples::getExample(99), std::invalid_argument);
}

// --- Example 3: Identity (m=4, all circles) ---

TEST(ThesisExamples, Example3HasCorrectConnectivity)
{
    auto preset = ThesisExamples::getExample(3);
    EXPECT_EQ(preset.target_domain->getConnectivity(), 4);
}

TEST(ThesisExamples, Example3HasCorrectInitialGuesses)
{
    auto preset = ThesisExamples::getExample(3);
    EXPECT_EQ(preset.initial_centers.size(), 3u);
    EXPECT_EQ(preset.initial_radii.size(), 3u);

    // First inner boundary initial guess: c=-0.4, rho=0.25
    EXPECT_NEAR(preset.initial_centers[0].real(), -0.4, 1e-12);
    EXPECT_NEAR(preset.initial_centers[0].imag(), 0.0, 1e-12);
    EXPECT_NEAR(preset.initial_radii[0], 0.25, 1e-12);
}

TEST(ThesisExamples, Example3HasCorrectConfig)
{
    auto preset = ThesisExamples::getExample(3);
    EXPECT_EQ(preset.config.N, 256);
    EXPECT_NEAR(preset.config.newton_tolerance, 1e-14, 1e-20);
    EXPECT_NEAR(preset.config.cgm_tolerance, 1e-15, 1e-20);
    EXPECT_EQ(preset.config.initial_guess_method, FornbergMCConfiguration::InitialGuessMethod::MANUAL);
}

TEST(ThesisExamples, Example3HasNameAndDescription)
{
    auto preset = ThesisExamples::getExample(3);
    EXPECT_FALSE(preset.name.empty());
    EXPECT_FALSE(preset.description.empty());
}

// --- Example 5: Ellipses (m=3) ---

TEST(ThesisExamples, Example5HasCorrectConnectivity)
{
    auto preset = ThesisExamples::getExample(5);
    EXPECT_EQ(preset.target_domain->getConnectivity(), 3);
}

TEST(ThesisExamples, Example5HasCorrectInitialGuesses)
{
    auto preset = ThesisExamples::getExample(5);
    EXPECT_EQ(preset.initial_centers.size(), 2u);
    EXPECT_EQ(preset.initial_radii.size(), 2u);
}

TEST(ThesisExamples, Example5HasCorrectN)
{
    auto preset = ThesisExamples::getExample(5);
    EXPECT_EQ(preset.config.N, 256);
}

// --- Example 2: Mixed (m=4, inverted ellipse outer) ---

TEST(ThesisExamples, Example2HasCorrectConnectivity)
{
    auto preset = ThesisExamples::getExample(2);
    EXPECT_EQ(preset.target_domain->getConnectivity(), 4);
}

TEST(ThesisExamples, Example2HasCorrectInitialGuesses)
{
    auto preset = ThesisExamples::getExample(2);
    EXPECT_EQ(preset.initial_centers.size(), 3u);
    EXPECT_EQ(preset.initial_radii.size(), 3u);
}

TEST(ThesisExamples, Example2HasCorrectN)
{
    auto preset = ThesisExamples::getExample(2);
    EXPECT_EQ(preset.config.N, 128);
}

// --- Example 4: High connectivity (m=7) ---

TEST(ThesisExamples, Example4HasCorrectConnectivity)
{
    auto preset = ThesisExamples::getExample(4);
    EXPECT_EQ(preset.target_domain->getConnectivity(), 7);
}

TEST(ThesisExamples, Example4HasCorrectInitialGuesses)
{
    auto preset = ThesisExamples::getExample(4);
    EXPECT_EQ(preset.initial_centers.size(), 6u);
    EXPECT_EQ(preset.initial_radii.size(), 6u);
}

TEST(ThesisExamples, Example4HasCorrectN)
{
    auto preset = ThesisExamples::getExample(4);
    EXPECT_EQ(preset.config.N, 128);
}

// --- Cross-cutting: initial guess count matches inner boundary count ---

TEST(ThesisExamples, AllPresetsHaveMatchingGuessAndBoundaryCount)
{
    for (int n : ThesisExamples::availableExamples())
    {
        auto preset = ThesisExamples::getExample(n);
        int inner_count = preset.target_domain->getConnectivity() - 1;
        EXPECT_EQ(static_cast<int>(preset.initial_centers.size()), inner_count)
            << "Example " << n << ": center count mismatch";
        EXPECT_EQ(static_cast<int>(preset.initial_radii.size()), inner_count)
            << "Example " << n << ": radii count mismatch";
    }
}

TEST(ThesisExamples, AllPresetBoundariesEvaluateToFinitePoints)
{
    for (int n : ThesisExamples::availableExamples())
    {
        auto preset = ThesisExamples::getExample(n);
        auto& boundaries = preset.target_domain->getBoundaries();
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            Complex z = boundaries[i]->evaluate(0.0);
            EXPECT_TRUE(std::isfinite(z.real()))
                << "Example " << n << ", boundary " << i << ": non-finite real part";
            EXPECT_TRUE(std::isfinite(z.imag()))
                << "Example " << n << ", boundary " << i << ": non-finite imag part";
        }
    }
}
