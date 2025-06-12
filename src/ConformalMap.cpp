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

#include "ConformalMap.h"
#include "ConformalMapMethod.h"

#include <stdexcept>

ConformalMap::ConformalMap(
    std::shared_ptr<Domain> source,
    std::shared_ptr<Domain> target,
    std::shared_ptr<ConformalMapMethod> method_impl
)
    : mp_source_domain(source)
    , mp_target_domain(target)
    , mp_method(method_impl)
    , m_mapping_type(determineMappingType(*source, *target))
{
    if (!source)
    {
        throw std::invalid_argument("Source domain cannot be null");
    }

    if (!target)
    {
        throw std::invalid_argument("Target domain cannot be null");
    }
}

void ConformalMap::setMethod(std::shared_ptr<ConformalMapMethod> method_impl)
{
    mp_method = method_impl;
}

void ConformalMap::compute(double target_accuracy)
{
    if (!mp_method)
    {
        throw std::runtime_error("No method set for computation");
    }

    // Validate domains before computation
    mp_method->validateDomains(*this);

    mp_method->compute(*this, target_accuracy);
}

Complex ConformalMap::map(const Complex& z) const
{
    if (!mp_method)
    {
        throw std::runtime_error("No method set for map evaluation");
    }

    return mp_method->map(z);
}

Complex ConformalMap::inverseMap(const Complex& w) const
{
    if (!mp_method)
    {
        throw std::runtime_error("No method set for inverse map evaluation");
    }

    return mp_method->inverseMap(w);
}

