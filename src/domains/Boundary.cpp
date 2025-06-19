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

#include "Boundary.h"

#include <stdexcept>

// Boundary implementations
Boundary::Boundary(std::vector<std::shared_ptr<BoundaryComponent>> boundaryComponents)
{
    for (const auto& component : boundaryComponents)
    {
        addComponent(component);
    }
}

Boundary::Boundary(std::shared_ptr<BoundaryComponent> component)
{
    addComponent(component);
}

void Boundary::addComponent(std::shared_ptr<BoundaryComponent> component)
{
    components.push_back(component);
}

const BoundaryComponent& Boundary::getComponent(size_t index) const
{
    if (index >= components.size())
    {
        throw std::out_of_range("Boundary component index out of range");
    }
    return *components[index];
}

Complex Boundary::evaluate(double t, size_t componentIndex) const
{
    if (components.empty())
    {
        throw std::runtime_error("Cannot evaluate empty boundary");
    }

    if (componentIndex >= components.size())
    {
        throw std::out_of_range("Component index out of range");
    }

    return components[componentIndex]->evaluate(t);
}

Complex Boundary::evaluateDerivative(double t, size_t componentIndex) const
{
    if (components.empty())
    {
        throw std::runtime_error("Cannot evaluate derivative of empty boundary");
    }

    if (componentIndex >= components.size())
    {
        throw std::out_of_range("Component index out of range");
    }

    return components[componentIndex]->evaluateDerivative(t);
}

std::vector<std::vector<Complex>> Boundary::sample(size_t numPoints) const
{
    if (components.empty() || numPoints == 0)
    {
        return {};
    }

    std::vector<std::vector<Complex>> samples;
    samples.reserve(components.size());

    for (size_t i = 0; i < components.size(); ++i)
    {
        samples.push_back(components[i]->sample(numPoints));
    }

    return samples;
}

double Boundary::findParameterization(const Complex& z, size_t componentIndex) const
{
    if (components.empty())
    {
        throw std::runtime_error("Cannot find parameterization on empty boundary");
    }

    if (componentIndex >= components.size())
    {
        throw std::out_of_range("Component index out of range");
    }

    return components[componentIndex]->findParameterization(z);
}


