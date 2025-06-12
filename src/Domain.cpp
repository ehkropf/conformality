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

#include "Domain.h"
#include <limits>
#include <algorithm>
#include <cmath>

using std::imag;
using std::real;

// Domain implementation
int Domain::calculateWindingNumber(const Complex& z, const std::vector<Complex>& samples, double tolerance) const
{
    if (samples.empty())
    {
        return 0;
    }

    int winding = 0;
    double minDistanceToBoundary = std::numeric_limits<double>::max();

    for (size_t i = 0; i < samples.size(); ++i)
    {
        Complex p1 = samples[i];
        Complex p2 = samples[(i + 1) % samples.size()];

        // Track minimum distance to boundary for tolerance checking
        double distToSegment = distanceToLineSegment(z, p1, p2);
        minDistanceToBoundary = std::min(minDistanceToBoundary, distToSegment);

        // Check if point is very close to boundary
        if (distToSegment < tolerance)
        {
            // Point is on boundary - return special value
            return std::numeric_limits<int>::max();
        }

        // Ray casting: check if horizontal ray from z to the right crosses this edge
        if ((imag(p1) <= imag(z) && imag(p2) > imag(z)) ||
            (imag(p2) <= imag(z) && imag(p1) > imag(z)))
        {
            // Avoid division by zero
            if (std::abs(imag(p2) - imag(p1)) < 1e-15)
            {
                continue;
            }

            // Calculate x-coordinate of intersection point
            double t = (imag(z) - imag(p1)) / (imag(p2) - imag(p1));
            double x_intersect = real(p1) + t * (real(p2) - real(p1));

            // Count intersection if it's to the right of the point
            if (x_intersect > real(z))
            {
                if (imag(p2) > imag(p1))
                {
                    winding++;
                }
                else
                {
                    winding--;
                }
            }
        }
    }

    return winding;
}

double Domain::distanceToLineSegment(const Complex& point, const Complex& segmentStart, const Complex& segmentEnd) const
{
    Complex v = segmentEnd - segmentStart;
    Complex w = point - segmentStart;

    // If segment has zero length, return distance to start point
    double segmentLengthSq = std::norm(v);
    if (segmentLengthSq < 1e-15)
    {
        return std::abs(point - segmentStart);
    }

    // Project point onto line segment
    double t = real(w * Complex(real(v), -imag(v))) / segmentLengthSq;
    t = std::max(0.0, std::min(1.0, t)); // Clamp to [0,1]

// FIXME:￼ Only need to scale by t no need for full complex
    Complex projection = segmentStart + v * Complex(t, 0.0);
    return std::abs(point - projection);
}

// SimplyConnectedDomain implementation
bool SimplyConnectedDomain::contains(const Complex& z) const
{
    // Increase sampling resolution for better accuracy (Key Change #3)
    const int sampleCount = 500; // Increased from 100
    auto samplesVec = boundary->sample(sampleCount);
    if (samplesVec.empty() || samplesVec[0].empty())
    {
        return false;
    }
    auto samples = samplesVec[0]; // Use first component for simply connected domain

// FIXME: should be somewhere else!!!!
    // Add proper boundary tolerance handling (Key Change #2)
    const double boundaryTolerance = 1e-12;

    // Fix the ray-casting intersection calculation and ensure correct winding number (Key Changes #1 and #4)
    int winding = calculateWindingNumber(z, samples, boundaryTolerance);

    // Handle boundary cases (Key Change #2)
    if (winding == std::numeric_limits<int>::max())
    {
        // Point is on or very close to boundary
        return !isExternal; // Inside for internal domains, outside for external domains
    }

    if (isExternal)
    {
        // For external domain, point is inside if outside the boundary (winding number = 0)
        return (winding == 0);
    }
    else
    {
        // For internal domain, point is inside if inside the boundary (winding number != 0)
        return (winding != 0);
    }
}

void SimplyConnectedDomain::transformBoundary(std::function<Complex(const Complex&)> transform)
{
    // Increase sampling resolution for transformed boundaries (Key Change #3)
    const int transformSampleCount = 2000; // Increased from 1000
    auto samplesVec = boundary->sample(transformSampleCount);
    std::vector<Complex> transformedSamples;
    auto samples = samplesVec[0]; // Use first component for simply connected domain
    transformedSamples.reserve(samples.size());

    for (const auto& point : samples)
    {
        transformedSamples.push_back(transform(point));
    }

    // Create a new boundary from transformed samples
    // Note: In a real implementation, we would need to handle
    // different boundary types more intelligently
    auto discreteComponent = std::make_shared<DiscreteBoundaryComponent>(transformedSamples);
    boundary = std::make_shared<Boundary>(discreteComponent);
}

// StarlikeDomain implementation
bool StarlikeDomain::contains(const Complex& z) const
{
    // For starlike domains, we can use the more efficient analytical containment test
    // with improved boundary tolerance handling (Key Change #2)
    const double boundaryTolerance = 1e-12;

    double angle = std::arg(z - center);
    double radius = std::abs(z - center);
    double boundaryRadius = radiusFunction(angle);

    // Handle boundary cases with tolerance
    if (std::abs(radius - boundaryRadius) < boundaryTolerance)
    {
        return !isExternal; // On boundary: inside for internal domains, outside for external domains
    }

    if (isExternal)
    {
        // For external domain, point is inside if it's outside the boundary
        return radius > boundaryRadius;
    }
    else
    {
        // For internal domain, point is inside if it's inside the boundary
        return radius < boundaryRadius;
    }
}

std::shared_ptr<Boundary> StarlikeDomain::createBoundary(
    const Complex& center,
    std::function<double(double)> radiusFunc
)
{
    // Create a parameterization from the radius function
    auto paramFunc = [center, radiusFunc](double t) -> Complex
    {
        double r = radiusFunc(t);
        return center + Complex(r * std::cos(t), r * std::sin(t));
    };

    // Compute the derivative of the parameterization
    auto derivFunc = [radiusFunc](double t) -> Complex
    {
        // For a radius function r(θ), the derivative of r(θ)e^(iθ) is:
        // r'(θ)e^(iθ) + ir(θ)e^(iθ)
        double r = radiusFunc(t);

        // Approximate r'(θ) with finite difference
        const double h = 1e-6;
        double rPrime = (radiusFunc(t + h) - radiusFunc(t - h)) / (2.0 * h);

        Complex tangent = Complex(rPrime * std::cos(t), rPrime * std::sin(t)) +
                                Complex(0.0, 1.0) * Complex(r * std::cos(t), r * std::sin(t));

        return tangent;
    };

    auto component = std::make_shared<AnalyticBoundaryComponent>(paramFunc, derivFunc);
    return std::make_shared<Boundary>(component);
}

// EllipticalDomain implementation

// CircularDomain implementation
bool CircularDomain::contains(const Complex& z) const
{
    const double boundaryTolerance = 1e-12;
    double dist = std::abs(z - getCenter());

    // Handle boundary cases with tolerance
    if (std::abs(dist - radius) < boundaryTolerance)
    {
        return !isExternalDomain(); // On boundary: inside for internal domains, outside for external domains
    }

    if (isExternalDomain())
    {
        return dist > radius;
    }
    else
    {
        return dist < radius;
    }
}

// PolygonalDomain implementation
void PolygonalDomain::setVertices(const std::vector<Complex>& newVertices)
{
    vertices = newVertices;
    boundary = createBoundary(vertices);
}

std::shared_ptr<Boundary> PolygonalDomain::createBoundary(
    const std::vector<Complex>& vertices
)
{
    // Create a piecewise linear parameterization from the vertices
    auto paramFunc = [vertices](double t) -> Complex
    {
        if (vertices.empty())
        {
            return Complex(0.0, 0.0);
        }

        double normalizedT = std::fmod(t, 2.0 * M_PI);
        if (normalizedT < 0) normalizedT += 2.0 * M_PI;

        double indexF = normalizedT * vertices.size() / (2.0 * M_PI);
        int index1 = static_cast<int>(std::floor(indexF)) % vertices.size();
        int index2 = (index1 + 1) % vertices.size();
        double frac = indexF - index1;

        return vertices[index1] * Complex(1.0 - frac) + vertices[index2] * Complex(frac);
    };

    // Compute the derivative of the parameterization
    auto derivFunc = [vertices](double t) -> Complex
    {
        if (vertices.empty())
        {
            return Complex(0.0, 0.0);
        }

        double normalizedT = std::fmod(t, 2.0 * M_PI);
        if (normalizedT < 0) normalizedT += 2.0 * M_PI;

        double indexF = normalizedT * vertices.size() / (2.0 * M_PI);
        int index1 = static_cast<int>(std::floor(indexF)) % vertices.size();
        int index2 = (index1 + 1) % vertices.size();

        // Derivative is constant along each edge
        return (vertices[index2] - vertices[index1]) * Complex(vertices.size() / (2.0 * M_PI));
    };

    auto component = std::make_shared<AnalyticBoundaryComponent>(paramFunc, derivFunc);
    return std::make_shared<Boundary>(component);
}

// MultiplyConnectedDomain implementation
void MultiplyConnectedDomain::addBoundary(std::shared_ptr<Boundary> boundary)
{
    boundaries.push_back(boundary);
    connectivity = boundaries.size();
}

bool MultiplyConnectedDomain::contains(const Complex& z) const
{
    if (boundaries.empty())
    {
        return false;
    }

    // For internal domains:
    // - Point is inside if it's inside the outer boundary (index 0)
    // - AND outside all inner boundaries (index > 0)
    //
    // For external domains:
    // - Point is inside if it's outside all boundaries

    bool insideOuterBoundary = false;
    const double boundaryTolerance = 1e-12;

    // Check against each boundary
    for (size_t i = 0; i < boundaries.size(); ++i)
    {
        // Increase sampling resolution for better accuracy (Key Change #3)
        auto samplesVec = boundaries[i]->sample(500);
        if (samplesVec.empty() || samplesVec[0].empty())
        {
            continue;
        }
        auto samples = samplesVec[0]; // Use first component of each boundary

        // Use improved winding number calculation (Key Changes #1 and #4)
        int winding = calculateWindingNumber(z, samples, boundaryTolerance);

        // Handle boundary cases (Key Change #2)
        if (winding == std::numeric_limits<int>::max())
        {
            // Point is on boundary
            if (i == 0)
            {
                // On outer boundary
                return !isExternal;
            }
            else
            {
                // On inner boundary
                return isExternal;
            }
        }

        bool insideThisBoundary = (winding != 0);

        if (i == 0)
        {
            // Outer boundary
            insideOuterBoundary = insideThisBoundary;

            // For external domains, being outside the outer boundary counts as "inside" the domain
            if (isExternal && !insideOuterBoundary)
            {
                return true;
            }

            // For internal domains, being outside the outer boundary means the point is outside the domain
            if (!isExternal && !insideOuterBoundary)
            {
                return false;
            }
        }
        else
        {
            // Inner boundary
            // For internal domains, being inside an inner boundary means the point is outside the domain
            if (!isExternal && insideThisBoundary)
            {
                return false;
            }

            // For external domains, being inside an inner boundary means the point is inside the domain
            if (isExternal && insideThisBoundary)
            {
                return true;
            }
        }
    }

    // If we get here for an internal domain, the point is inside the outer boundary
    // and outside all inner boundaries, so it's inside the domain
    if (!isExternal)
    {
        return true;
    }

    // If we get here for an external domain, the point is inside the outer boundary
    // and outside all inner boundaries, so it's outside the domain
    return false;
}

void MultiplyConnectedDomain::transformBoundary(std::function<Complex(const Complex&)> transform)
{
    // Transform each boundary with increased sampling resolution (Key Change #3)
    const int transformSampleCount = 2000;

    for (auto& boundary : boundaries)
    {
        auto samplesVec = boundary->sample(transformSampleCount);
        std::vector<Complex> transformedSamples;
        auto samples = samplesVec[0]; // Use first component of each boundary
        transformedSamples.reserve(samples.size());

        for (const auto& point : samples)
        {
            transformedSamples.push_back(transform(point));
        }

        // Create a new boundary from transformed samples
        auto discreteComponent = std::make_shared<DiscreteBoundaryComponent>(transformedSamples);
        boundary = std::make_shared<Boundary>(discreteComponent);
    }
}
