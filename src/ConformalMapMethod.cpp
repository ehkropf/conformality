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

#include "ConformalMapMethod.h"
#include "ConformalMap.h"
#include "Domain.h"
#include <stdexcept>

ConformalMapMethod::ConformalMapMethod()
    : m_achieved_accuracy(0.0)
    , m_max_iterations(1000)
    , m_iteration_count(0)
{
}

ConformalMapMethod::ConformalMapMethod(int max_iter)
    : m_achieved_accuracy(0.0)
    , m_max_iterations(max_iter)
    , m_iteration_count(0)
{
    if (max_iter <= 0)
    {
        throw std::invalid_argument("Maximum iterations must be positive");
    }
}

void ConformalMapMethod::setMaxIterations(int max)
{
    if (max <= 0)
    {
        throw std::invalid_argument("Maximum iterations must be positive");
    }

    m_max_iterations = max;
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

void ConformalMapMethod::validateDomain(
    std::shared_ptr<Domain> domain,
    int expected_connectivity
) const
{
    validateDomainCompatibility(domain, expected_connectivity);
    validateDomainGeometry(domain);
}

void ConformalMapMethod::validateDomains(const ConformalMap& map_instance) const
{
    const auto source_domain = map_instance.getSourceDomain();
    const auto target_domain = map_instance.getTargetDomain();

    if (!source_domain || !target_domain)
    {
        throw std::invalid_argument("Both source and target domains must be non-null");
    }

    int source_connectivity = source_domain->getConnectivity();
    int target_connectivity = target_domain->getConnectivity();

    if (source_connectivity != target_connectivity)
    {
        throw std::invalid_argument(
            "Source and target domains must have the same connectivity. "
            "Source: " + std::to_string(source_connectivity) +
            ", Target: " + std::to_string(target_connectivity)
        );
    }

    // Validate individual domain geometry
    validateDomainGeometry(source_domain);
    validateDomainGeometry(target_domain);
}
