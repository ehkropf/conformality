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

#pragma once

#include "../domains/Domain.h"
#include "../methods/FornbergMCConfiguration.h"

#include <memory>
#include <string>
#include <vector>

namespace conformality::examples
{

/// Preset configuration for a thesis example, including domain, initial guesses, and solver config.
struct ThesisExamplePreset
{
    std::string name;
    std::string description;
    std::shared_ptr<MultiplyConnectedDomain> target_domain;
    std::vector<Complex> initial_centers;
    std::vector<double> initial_radii;
    FornbergMCConfiguration config;
};

/// Provides hardcoded thesis example presets (from MATLAB th_gen_ex*.m files).
class ThesisExamples
{
public:
    static ThesisExamplePreset getExample(int exampleNumber);
    static std::vector<int> availableExamples();

private:
    static ThesisExamplePreset makeExample1();
    static ThesisExamplePreset makeExample2();
    static ThesisExamplePreset makeExample3();
    static ThesisExamplePreset makeExample4();
    static ThesisExamplePreset makeExample5();
    static FornbergMCConfiguration makeConfig(int N);
};

} // namespace conformality::examples
