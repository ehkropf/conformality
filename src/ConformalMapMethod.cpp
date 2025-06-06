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
#include "Domain.h"
#include <stdexcept>

ConformalMapMethod::ConformalMapMethod()
    : accuracy(0.0)
    , max_iterations(1000)
    , iteration_count(0)
{
}

void ConformalMapMethod::setMaxIterations(int max)
{
    if (max <= 0)
    {
        throw std::invalid_argument("Maximum iterations must be positive");
    }

    max_iterations = max;
}

void ConformalMapMethod::validateDomainCompatibility(
    std::shared_ptr<Domain> domain,
    int expected_connectivity
) const
{
    if (!domain)
    {
        throw std::invalid_argument("Domain cannot be null");
    }

    if (domain->getConnectivity() != expected_connectivity)
    {
        throw std::invalid_argument(
            "Domain connectivity mismatch. Expected: " + std::to_string(expected_connectivity) +
            ", Actual: " + std::to_string(domain->getConnectivity())
        );
    }
}

void ConformalMapMethod::validateDomainGeometry(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("Domain cannot be null");
    }

    // This method should be overridden by specific methods to check
    // for domain-specific geometric properties (e.g., starlike, polygonal)
    // Base implementation just validates non-null domain
    // Derived classes can use dynamic_cast to check specific domain types
}
