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
#include "../src/core/MethodInfo.h"
#include "../src/methods/ConformalMapMethod.h"
#include "../src/methods/FornbergMC.h"
#include "../src/methods/FornbergMCConfiguration.h"
#include "../src/methods/PMatrixBuilder.h"
#include "../src/numerics/CGSolver.h"

using conformality::formatMethodInfoValue;
using conformality::MethodInfo;
using conformality::MethodInfoValue;

// ---------- Base class default ----------

TEST(MethodInfoBaseClassTest, ReturnsEmpty)
{
    // Verify default-constructed MethodInfo has empty fields
    MethodInfo empty{};
    EXPECT_TRUE(empty.name.empty());
    EXPECT_TRUE(empty.parameters.empty());
    EXPECT_TRUE(empty.results.empty());
}

// ---------- FornbergMC parameters ----------

TEST(MethodInfoFornbergMCTest, ParametersMatchConfig)
{
    FornbergMCConfiguration config;
    config.N = 128;
    config.newton_tolerance = 1e-12;
    config.max_newton_iterations = 50;
    FornbergMC fm(config);

    auto info = fm.getMethodInfo();
    EXPECT_EQ(info.name, "FornbergMC");
    ASSERT_EQ(info.parameters.size(), 5u);

    // Connectivity (default 0 before compute)
    EXPECT_EQ(info.parameters[0].label, "Connectivity");
    EXPECT_EQ(std::get<int>(info.parameters[0].value), 0);

    // N
    EXPECT_EQ(info.parameters[1].label, "N");
    EXPECT_EQ(std::get<int>(info.parameters[1].value), 128);

    // Newton tolerance
    EXPECT_EQ(info.parameters[2].label, "Newton tolerance");
    EXPECT_DOUBLE_EQ(std::get<double>(info.parameters[2].value), 1e-12);

    // Max Newton iterations
    EXPECT_EQ(info.parameters[3].label, "Max Newton iterations");
    EXPECT_EQ(std::get<int>(info.parameters[3].value), 50);

    // Annulus case (default false before compute)
    EXPECT_EQ(info.parameters[4].label, "Annulus case");
    EXPECT_FALSE(std::get<bool>(info.parameters[4].value));
}

// ---------- FornbergMC results after state setup ----------

TEST(MethodInfoFornbergMCTest, ResultsAfterStateSetup)
{
    FornbergMC fm;

    // Set internal state via FRIEND_TEST access
    fm.m_is_converged = true;
    fm.m_current_residual = 1.5e-10;
    fm.m_residual_history = {1.0, 0.5, 0.1, 1.5e-10};

    auto info = fm.getMethodInfo();
    ASSERT_EQ(info.results.size(), 3u);

    EXPECT_EQ(info.results[0].label, "Iterations");
    EXPECT_EQ(std::get<int>(info.results[0].value), 4);

    EXPECT_EQ(info.results[1].label, "Residual");
    EXPECT_DOUBLE_EQ(std::get<double>(info.results[1].value), 1.5e-10);

    EXPECT_EQ(info.results[2].label, "Converged");
    EXPECT_TRUE(std::get<bool>(info.results[2].value));
}

// ---------- formatMethodInfoValue ----------

TEST(MethodInfoFormatTest, FormatsInt)
{
    MethodInfoValue v = 42;
    EXPECT_EQ(formatMethodInfoValue(v), "42");
}

TEST(MethodInfoFormatTest, FormatsDouble)
{
    MethodInfoValue v = 1.5e-10;
    auto s = formatMethodInfoValue(v);
    EXPECT_NE(s.find("1.50e"), std::string::npos);
}

TEST(MethodInfoFormatTest, FormatsBool)
{
    EXPECT_EQ(formatMethodInfoValue(MethodInfoValue{true}), "Yes");
    EXPECT_EQ(formatMethodInfoValue(MethodInfoValue{false}), "No");
}

TEST(MethodInfoFormatTest, FormatsString)
{
    MethodInfoValue v = std::string("hello");
    EXPECT_EQ(formatMethodInfoValue(v), "hello");
}

// ---------- Type safety ----------

TEST(MethodInfoTypeSafetyTest, VariantHoldsExpectedTypes)
{
    FornbergMC fm;
    fm.m_is_converged = false;
    fm.m_current_residual = 2.0e-5;
    fm.m_residual_history = {1.0, 2.0e-5};

    auto info = fm.getMethodInfo();

    // Verify type tags
    EXPECT_TRUE(std::holds_alternative<int>(info.results[0].value));       // Iterations
    EXPECT_TRUE(std::holds_alternative<double>(info.results[1].value));    // Residual
    EXPECT_TRUE(std::holds_alternative<bool>(info.results[2].value));      // Converged

    // Verify std::get works without throwing
    EXPECT_NO_THROW(std::get<int>(info.results[0].value));
    EXPECT_NO_THROW(std::get<double>(info.results[1].value));
    EXPECT_NO_THROW(std::get<bool>(info.results[2].value));
}
