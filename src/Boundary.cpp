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

#include <cmath>
#include <limits>
#include <stdexcept>

// Boundary implementations
void Boundary::addComponent(std::shared_ptr<BoundaryComponent> component)
{
    component->setIndex(components.size());
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


// SimpleBoundary implementations
SimpleBoundary::SimpleBoundary(std::shared_ptr<BoundaryComponent> component)
{
    addComponent(component);
}

// CompositeBoundary implementations
CompositeBoundary::CompositeBoundary(std::vector<std::shared_ptr<BoundaryComponent>> boundaryComponents)
{
    for (auto& component : boundaryComponents)
    {
        addComponent(component);
    }
    updateParameterRanges();
}

void CompositeBoundary::addComponent(std::shared_ptr<BoundaryComponent> component)
{
    Boundary::addComponent(component);
    updateParameterRanges();
}

Complex CompositeBoundary::evaluate(double t) const
{
    double normalizedT = std::fmod(t, 2.0 * M_PI);
    if (normalizedT < 0) normalizedT += 2.0 * M_PI;

    // Find which component the parameter corresponds to
    size_t componentIdx = 0;
    double componentT = normalizedT;

    for (size_t i = 0; i < parameterRanges.size(); ++i)
    {
        if (normalizedT < parameterRanges[i])
        {
            componentIdx = i;
            break;
        }
        componentT -= parameterRanges[i];
    }

    // Scale parameter to component's range [0, 2π)
    double scaledT = componentT * 2.0 * M_PI / parameterRanges[componentIdx];
    return components[componentIdx]->evaluate(scaledT);
}

Complex CompositeBoundary::evaluateDerivative(double t) const
{
    double normalizedT = std::fmod(t, 2.0 * M_PI);
    if (normalizedT < 0) normalizedT += 2.0 * M_PI;

    // Find which component the parameter corresponds to
    size_t componentIdx = 0;
    double componentT = normalizedT;

    for (size_t i = 0; i < parameterRanges.size(); ++i)
    {
        if (normalizedT < parameterRanges[i])
        {
            componentIdx = i;
            break;
        }
        componentT -= parameterRanges[i];
    }

    // Scale parameter to component's range [0, 2π)
    double scaledT = componentT * 2.0 * M_PI / parameterRanges[componentIdx];
    return components[componentIdx]->evaluateDerivative(scaledT);
}

std::vector<Complex> CompositeBoundary::sample(int numPoints) const
{
    if (components.empty() || numPoints <= 0)
    {
        return {};
    }

    std::vector<Complex> samples;
    samples.reserve(numPoints);

    // Distribute points proportionally to parameter ranges
    for (size_t i = 0; i < components.size(); ++i)
    {
        int componentPoints = static_cast<int>(numPoints * parameterRanges[i] / (2.0 * M_PI));
        if (componentPoints <= 0)
        {
            componentPoints = 1;  // Ensure at least one point per component
        }

        auto componentSamples = components[i]->sample(componentPoints);
        samples.insert(samples.end(), componentSamples.begin(), componentSamples.end());
    }

    // Adjust to exactly numPoints if needed
    if (samples.size() > static_cast<size_t>(numPoints))
    {
        samples.resize(numPoints);
    }
    else if (samples.size() < static_cast<size_t>(numPoints))
    {
        // Add more points to the first component if needed
        int remaining = numPoints - samples.size();
        auto extraSamples = components[0]->sample(remaining);
        samples.insert(samples.end(), extraSamples.begin(), extraSamples.end());
    }

    return samples;
}

double CompositeBoundary::findParameterization(const Complex& z) const
{
    if (components.empty())
    {
        return 0.0;
    }

    // Find closest component and its parameter
    double minDist = std::numeric_limits<double>::max();
    double bestParam = 0.0;

    for (size_t i = 0; i < components.size(); ++i)
    {
        double t = components[i]->findParameterization(z);
        Complex point = components[i]->evaluate(t);
        double dist = std::norm(z - point);

        if (dist < minDist)
        {
            minDist = dist;
            // Map local parameter to global parameter
            double globalParam = 0.0;
            for (size_t j = 0; j < i; ++j)
            {
                globalParam += parameterRanges[j];
            }
            globalParam += t * parameterRanges[i] / (2.0 * M_PI);
            bestParam = globalParam;
        }
    }

    return bestParam;
}

void CompositeBoundary::updateParameterRanges()
{
    parameterRanges.resize(components.size());

    // Simple approach: divide parameter space equally
    double paramPerComponent = 2.0 * M_PI / components.size();

    for (size_t i = 0; i < components.size(); ++i)
    {
        parameterRanges[i] = paramPerComponent;
    }

    // Convert to cumulative ranges
    for (size_t i = 1; i < parameterRanges.size(); ++i)
    {
        parameterRanges[i] += parameterRanges[i-1];
    }
}