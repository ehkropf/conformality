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

#include "ConformalMap.h"
#include "ConformalMapMethod.h"
#include <stdexcept>

ConformalMap::ConformalMap(
    std::shared_ptr<Domain> source,
    std::shared_ptr<Domain> target,
    bool external
)
    : source_domain(source)
    , target_domain(target)
    , is_external(external)
    , method(nullptr)
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

std::shared_ptr<Domain> ConformalMap::getSourceDomain() const
{
    return source_domain;
}

std::shared_ptr<Domain> ConformalMap::getTargetDomain() const
{
    return target_domain;
}

bool ConformalMap::isExternalMap() const
{
    return is_external;
}

void ConformalMap::setExternal(bool external)
{
    is_external = external;
}

void ConformalMap::setMethod(std::shared_ptr<ConformalMapMethod> method_impl)
{
    method = method_impl;
}

std::shared_ptr<ConformalMapMethod> ConformalMap::getMethod() const
{
    return method;
}

void ConformalMap::compute(double target_accuracy)
{
    if (!method)
    {
        throw std::runtime_error("No method set for computation");
    }

    method->compute(*this, target_accuracy);
}
