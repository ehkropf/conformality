/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
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

#include "ConformalMapMethod.h"
#include <stdexcept>
#include <typeinfo>

ConformalMapMethod::ConformalMapMethod()
    : accuracy(0.0)
    , max_iterations(1000)
    , iteration_count(0)
{
}

double ConformalMapMethod::getAccuracy() const
{
    return accuracy;
}

int ConformalMapMethod::getIterationCount() const
{
    return iteration_count;
}

void ConformalMapMethod::setMaxIterations(int max)
{
    if (max <= 0)
    {
        throw std::invalid_argument("Maximum iterations must be positive");
    }

    max_iterations = max;
}

void ConformalMapMethod::validateMapType(
    ConformalMap& map_instance,
    const std::string& expected_type
) const
{
    const std::string& actual_type = typeid(map_instance).name();

    // Check if the actual type contains the expected type
    // This is a simple check that works in many cases, but might need
    // to be more sophisticated depending on the exact typeid implementation
    if (actual_type.find(expected_type) == std::string::npos)
    {
        throw std::invalid_argument(
            "Map type mismatch. Expected: " + expected_type +
            ", Actual: " + actual_type
        );
    }
}
