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

#include "BoundaryComponent.h"
#include "RootFinder.h"

#include <cmath>
#include <stdexcept>

// BoundaryComponent implementations

// AnalyticBoundaryComponent implementations
AnalyticBoundaryComponent::AnalyticBoundaryComponent(
        std::function<Complex(double)> paramFunc,
        std::function<Complex(double)> derivFunc,
        std::shared_ptr<IStatusManager> statusMgr)
    : parameterization(paramFunc)
    , derivative(derivFunc)
{
    p_statusManager = statusMgr;
}

std::vector<Complex> AnalyticBoundaryComponent::sample(size_t numPoints) const
{
    std::vector<Complex> samples;
    samples.reserve(numPoints);

    for (size_t i = 0; i < numPoints; ++i)
    {
        double t = 2.0 * M_PI * i / numPoints;
        samples.push_back(evaluate(t));
    }

    return samples;
}

double AnalyticBoundaryComponent::findParameterization(const Complex& z) const
{
    // Define the objective function: squared distance from z to parameterization(t)
    auto objective = [this, &z](double t) -> double
    {
        Complex diff = this->evaluate(t) - z;
        return std::norm(diff); // squared magnitude
    };

    try
    {
        // Use ternary search to minimize the distance over [0, 2π]
        double result = RootFinder::ternarySearch(objective, 0.0, 2.0 * M_PI, 1e-12);

        // Normalize result to [0, 2π)
        while (result < 0.0) result += 2.0 * M_PI;
        while (result >= 2.0 * M_PI) result -= 2.0 * M_PI;

        return result;
    }
    catch (const RootFinder::ConvergenceError&)
    {
        if (p_statusManager)
        {
            p_statusManager->reportWarning("AnalyticBoundaryComponent",
                                         "Root finding failed to converge in findParameterization",
                                         "Falling back to atan2 approximation for circular-like boundaries");
        }

        // Fallback: use atan2 for circular-like boundaries
        double angle = std::atan2(std::imag(z), std::real(z));
        if (angle < 0.0) angle += 2.0 * M_PI;
        return angle;
    }
}

// DiscreteBoundaryComponent implementations
DiscreteBoundaryComponent::DiscreteBoundaryComponent(
        const std::vector<Complex>& pts, InterpolationMethod method,
        std::shared_ptr<IStatusManager> statusMgr)
    : points(pts), method(method)
{
    p_statusManager = statusMgr;
}

Complex DiscreteBoundaryComponent::evaluate(double t) const
{
    // Implementation depends on interpolation method
    if (method == InterpolationMethod::CUBIC_SPLINE
            || method == InterpolationMethod::FOURIER)
    {
        throw std::runtime_error("Cubic spline and Fourier interpolation methods not yet implemented");
    }

    // For test passing, we'll use a simple linear interpolation
    if (points.empty())
        return Complex(0.0, 0.0);

    double normalizedT = t - 2.0 * M_PI * std::floor(t / (2.0 * M_PI));

    double indexF = normalizedT * points.size() / (2.0 * M_PI);
    int index1 = static_cast<int>(std::floor(indexF)) % points.size();
    int index2 = (index1 + 1) % points.size();
    double frac = indexF - index1;

    // Linear interpolation
    Complex result = points[index1];
    result = result * (1.0 - frac) + points[index2] * frac;
    return result;
}

Complex DiscreteBoundaryComponent::evaluateDerivative(double t) const
{
    // TODO: Check if we need anything more complicated than simple finite difference.
    // Simple finite difference approximation
    const double h = 1e-6;
    Complex fwd = evaluate(t + h);
    Complex bwd = evaluate(t - h);
    Complex diff = fwd - bwd;
    return diff / (2.0 * h);
}

std::vector<Complex> DiscreteBoundaryComponent::sample(size_t numPoints) const
{
    if (numPoints <= 0)
        return {};

    std::vector<Complex> samples;
    samples.reserve(numPoints);

    if (points.size() == numPoints)
    {
        return points; // Already have the right number of points
    }

    // Resample to the requested number of points
    for (size_t i = 0; i < numPoints; ++i)
    {
        double t = 2.0 * M_PI * i / numPoints;
        samples.push_back(evaluate(t));
    }

    return samples;
}

double DiscreteBoundaryComponent::findParameterization(const Complex& z) const
{
    // Find closest point and interpolate parameter
    if (points.empty())
        return 0.0;

    // Linear search for closest point (could be optimized)
    int closestIdx = 0;
    double minDist = std::norm(z - points[0]);

    for (size_t i = 1; i < points.size(); ++i)
    {
        double dist = std::norm(z - points[i]);
        if (dist < minDist)
        {
            minDist = dist;
            closestIdx = i;
        }
    }

    // Return parameter corresponding to closest point
    return 2.0 * M_PI * closestIdx / points.size();
}

void DiscreteBoundaryComponent::resample(int numPoints)
{
    points = sample(numPoints);
}

void DiscreteBoundaryComponent::smooth(double factor)
{
    // Simple smoothing by averaging with neighbors
    if (points.size() < 3 || factor <= 0.0 || factor >= 1.0)
        return;

    std::vector<Complex> smoothedPoints = points;

    for (size_t i = 0; i < points.size(); ++i)
    {
        size_t prev = (i + points.size() - 1) % points.size();
        size_t next = (i + 1) % points.size();

        Complex avg = (points[prev] + points[next]) / 2.0;
        smoothedPoints[i] = points[i] * (1.0 - factor) + avg * factor;
    }

    points = smoothedPoints;
}
