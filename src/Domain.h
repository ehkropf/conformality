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

#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Boundary.h"
#include <functional>
#include <memory>
#include <cmath>

/**
 * @brief Abstract base class for domains
 */
class Domain
{
protected:
    bool isExternal;
    int connectivity;

    /**
     * @brief Improved ray-casting algorithm for point-in-polygon testing
     * @param z Point to test
     * @param samples Boundary sample points
     * @param tolerance Tolerance for boundary proximity
     * @return Winding number
     */
    int calculateWindingNumber(const ComplexDouble& z, const std::vector<ComplexDouble>& samples, double tolerance = 1e-12) const
    {
        if (samples.empty())
        {
            return 0;
        }

        int winding = 0;
        double minDistanceToBoundary = std::numeric_limits<double>::max();

        for (size_t i = 0; i < samples.size(); ++i)
        {
            ComplexDouble p1 = samples[i];
            ComplexDouble p2 = samples[(i + 1) % samples.size()];

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
            if ((p1.imag() <= z.imag() && p2.imag() > z.imag()) ||
                (p2.imag() <= z.imag() && p1.imag() > z.imag()))
            {
                // Avoid division by zero
                if (std::abs(p2.imag() - p1.imag()) < 1e-15)
                {
                    continue;
                }

                // Calculate x-coordinate of intersection point
                double t = (z.imag() - p1.imag()) / (p2.imag() - p1.imag());
                double x_intersect = p1.real() + t * (p2.real() - p1.real());

                // Count intersection if it's to the right of the point
                if (x_intersect > z.real())
                {
                    if (p2.imag() > p1.imag())
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

    /**
     * @brief Calculate distance from point to line segment
     * @param point Point to measure distance from
     * @param segmentStart Start of line segment
     * @param segmentEnd End of line segment
     * @return Distance to closest point on segment
     */
    double distanceToLineSegment(const ComplexDouble& point, const ComplexDouble& segmentStart, const ComplexDouble& segmentEnd) const
    {
        ComplexDouble v = segmentEnd - segmentStart;
        ComplexDouble w = point - segmentStart;

        // If segment has zero length, return distance to start point
        double segmentLengthSq = std::norm(v.getValue());
        if (segmentLengthSq < 1e-15)
        {
            return std::abs((point - segmentStart).getValue());
        }

        // Project point onto line segment
        double t = std::real((w * ComplexDouble(v.real(), -v.imag())).getValue()) / segmentLengthSq;
        t = std::max(0.0, std::min(1.0, t)); // Clamp to [0,1]

        ComplexDouble projection = segmentStart + v * ComplexDouble(t, 0.0);
        return std::abs((point - projection).getValue());
    }

public:
    /**
     * @brief Construct a new Domain
     * @param isExternalDomain Whether this is an external domain
     * @param domainConnectivity Connectivity of the domain
     */
    Domain(bool isExternalDomain = false, int domainConnectivity = 1)
        : isExternal(isExternalDomain)
        , connectivity(domainConnectivity)
    {}

    virtual ~Domain() = default;

    /**
     * @brief Check if this is an external domain
     * @return true if external, false if internal
     */
    bool isExternalDomain() const
    {
        return isExternal;
    }

    /**
     * @brief Get the connectivity of the domain
     * @return Connectivity (number of boundary components)
     */
    int getConnectivity() const
    {
        return connectivity;
    }

    /**
     * @brief Check if a point is contained in the domain
     * @param z ComplexDouble point to check
     * @return true if inside, false otherwise
     */
    virtual bool contains(const ComplexDouble& z) const = 0;

    /**
     * @brief Transform the boundary using a function
     * @param transform Function mapping complex points to complex points
     */
    virtual void transformBoundary(std::function<ComplexDouble(const ComplexDouble&)> transform) = 0;
};

/**
 * @brief Simply connected domain (connectivity = 1)
 */
class SimplyConnectedDomain : public Domain
{
protected:
    std::shared_ptr<Boundary> boundary;

public:
    /**
     * @brief Construct a new Simply Connected Domain
     * @param domainBoundary Boundary of the domain
     * @param isExternalDomain Whether this is an external domain
     */
    SimplyConnectedDomain(std::shared_ptr<Boundary> domainBoundary, bool isExternalDomain = false)
        : Domain(isExternalDomain, 1)
        , boundary(domainBoundary)
    {}

    Boundary& getBoundary()
    {
        return *boundary;
    }

    const Boundary& getBoundary() const
    {
        return *boundary;
    }

    bool contains(const ComplexDouble& z) const override
    {
        // Increase sampling resolution for better accuracy (Key Change #3)
        const int sampleCount = 500; // Increased from 100
        auto samples = boundary->sample(sampleCount);
        if (samples.empty())
        {
            return false;
        }

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

    void transformBoundary(std::function<ComplexDouble(const ComplexDouble&)> transform) override
    {
        // Increase sampling resolution for transformed boundaries (Key Change #3)
        const int transformSampleCount = 2000; // Increased from 1000
        auto samples = boundary->sample(transformSampleCount);
        std::vector<ComplexDouble> transformedSamples;
        transformedSamples.reserve(samples.size());

        for (const auto& point : samples)
        {
            transformedSamples.push_back(transform(point));
        }

        // Create a new boundary from transformed samples
        // Note: In a real implementation, we would need to handle
        // different boundary types more intelligently
        auto discreteComponent = std::make_shared<DiscreteBoundaryComponent>(transformedSamples);
        boundary = std::make_shared<SimpleBoundary>(discreteComponent);
    }
};

/**
 * @brief Starlike domain (simply connected domain that is star-shaped with respect to a center point)
 */
class StarlikeDomain : public SimplyConnectedDomain
{
protected:
    ComplexDouble center;
    std::function<double(double)> radiusFunction;

public:
    /**
     * @brief Construct a new Starlike Domain
     * @param domainCenter Center point of the domain
     * @param radiusFunc Function mapping angle to radius
     * @param isExternalDomain Whether this is an external domain
     */
    StarlikeDomain(
        const ComplexDouble& domainCenter,
        std::function<double(double)> radiusFunc,
        bool isExternalDomain = false
    )
        : SimplyConnectedDomain(createBoundary(domainCenter, radiusFunc), isExternalDomain)
        , center(domainCenter)
        , radiusFunction(radiusFunc)
    {}

    /**
     * @brief Get the center of the domain
     * @return Center point
     */
    ComplexDouble getCenter() const
    {
        return center;
    }

    /**
     * @brief Get the radius at a specific angle
     * @param angle Angle in radians
     * @return Radius at the given angle
     */
    double getRadius(double angle) const
    {
        return radiusFunction(angle);
    }

    bool contains(const ComplexDouble& z) const override
    {
        // For starlike domains, we can use the more efficient analytical containment test
        // with improved boundary tolerance handling (Key Change #2)
        const double boundaryTolerance = 1e-12;

        double angle = std::arg(z.getValue() - center.getValue());
        double radius = std::abs(z.getValue() - center.getValue());
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

private:
    static std::shared_ptr<Boundary> createBoundary(
        const ComplexDouble& center,
        std::function<double(double)> radiusFunc
    )
    {
        // Create a parameterization from the radius function
        auto paramFunc = [center, radiusFunc](double t) -> ComplexDouble
        {
            double r = radiusFunc(t);
            return center + ComplexDouble(r * std::cos(t), r * std::sin(t));
        };

        // Compute the derivative of the parameterization
        auto derivFunc = [radiusFunc](double t) -> ComplexDouble
        {
            // For a radius function r(θ), the derivative of r(θ)e^(iθ) is:
            // r'(θ)e^(iθ) + ir(θ)e^(iθ)
            double r = radiusFunc(t);

            // Approximate r'(θ) with finite difference
            const double h = 1e-6;
            double rPrime = (radiusFunc(t + h) - radiusFunc(t - h)) / (2.0 * h);

            ComplexDouble tangent = ComplexDouble(rPrime * std::cos(t), rPrime * std::sin(t)) +
                                    ComplexDouble(0.0, 1.0) * ComplexDouble(r * std::cos(t), r * std::sin(t));

            return tangent;
        };

        auto component = std::make_shared<AnalyticBoundaryComponent>(paramFunc, derivFunc);
        return std::make_shared<SimpleBoundary>(component);
    }
};

/**
 * @brief Elliptical domain
 */
class EllipticalDomain : public StarlikeDomain
{
private:
    double a;  // Semi-major axis
    double b;  // Semi-minor axis
    double rotation;  // Rotation angle in radians

public:
    /**
     * @brief Construct a new Elliptical Domain
     * @param semiMajorAxis Semi-major axis
     * @param semiMinorAxis Semi-minor axis
     * @param rotationAngle Rotation angle in radians
     * @param domainCenter Center of the ellipse
     * @param isExternalDomain Whether this is an external domain
     */
    EllipticalDomain(double semiMajorAxis, double semiMinorAxis, double rotationAngle = 0.0,
                     const ComplexDouble& domainCenter = ComplexDouble(0.0, 0.0),
                     bool isExternalDomain = false)
        : StarlikeDomain(domainCenter,
                         [semiMajorAxis, semiMinorAxis, rotationAngle, this](double angle) -> double
                         {
                             // Adjust angle by rotation.
                             double adjustedAngle = angle - rotationAngle;

                             const double theta = std::atan((semiMajorAxis*semiMinorAxis)*std::tan(adjustedAngle));
                             return std::abs(std::complex<double>(a*std::cos(theta), b*std::sin(theta)));

                         },
            isExternalDomain)
        , a(semiMajorAxis)
        , b(semiMinorAxis)
        , rotation(rotationAngle)
    {}

    /**
     * @brief Get the eccentricity of the ellipse
     * @return Eccentricity
     */
    double getEccentricity() const
    {
        if (a >= b)
        {
            return std::sqrt(1.0 - (b * b) / (a * a));
        }
        else
        {
            return std::sqrt(1.0 - (a * a) / (b * b));
        }
    }

    /**
     * @brief Get the semi-major axis
     * @return Semi-major axis
     */
    double getSemiMajorAxis() const
    {
        return a;
    }

    /**
     * @brief Get the semi-minor axis
     * @return Semi-minor axis
     */
    double getSemiMinorAxis() const
    {
        return b;
    }

    /**
     * @brief Get the rotation angle
     * @return Rotation angle in radians
     */
    double getRotation() const
    {
        return rotation;
    }
};

/**
 * @brief Circular domain (special case of elliptical domain with equal axes)
 */
class CircularDomain : public EllipticalDomain
{
private:
    double radius;

public:
    /**
     * @brief Construct a new Circular Domain
     * @param domainCenter Center of the circle
     * @param circleRadius Radius of the circle
     * @param isExternalDomain Whether this is an external domain
     */
    CircularDomain(
        const ComplexDouble& domainCenter,
        double circleRadius,
        bool isExternalDomain = false
    )
        : EllipticalDomain(
            circleRadius,
            circleRadius,
            0.0,
            domainCenter,
            isExternalDomain
        )
        , radius(circleRadius)
    {}

    /**
     * @brief Get the radius of the circle
     * @return Radius
     */
    double getRadius() const
    {
        return radius;
    }

    // Override contains for efficiency with improved boundary tolerance (Key Change #2)
    bool contains(const ComplexDouble& z) const override
    {
        const double boundaryTolerance = 1e-12;
        double dist = std::abs(z.getValue() - getCenter().getValue());

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
};

/**
 * @brief Polygonal domain (domain bounded by a polygon)
 */
class PolygonalDomain : public SimplyConnectedDomain
{
private:
    std::vector<ComplexDouble> vertices;

public:
    /**
     * @brief Construct a new Polygonal Domain
     * @param domainVertices Vertices of the polygon
     * @param isExternalDomain Whether this is an external domain
     */
    PolygonalDomain(
        const std::vector<ComplexDouble>& domainVertices,
        bool isExternalDomain = false
    )
        : SimplyConnectedDomain(createBoundary(domainVertices), isExternalDomain)
        , vertices(domainVertices)
    {}

    /**
     * @brief Get the vertices of the polygon
     * @return Vector of vertices
     */
    const std::vector<ComplexDouble>& getVertices() const
    {
        return vertices;
    }

    /**
     * @brief Set the vertices of the polygon
     * @param newVertices New vertices
     */
    void setVertices(const std::vector<ComplexDouble>& newVertices)
    {
        vertices = newVertices;
        boundary = createBoundary(vertices);
    }

private:
    static std::shared_ptr<Boundary> createBoundary(
        const std::vector<ComplexDouble>& vertices
    )
    {
        // Create a piecewise linear parameterization from the vertices
        auto paramFunc = [vertices](double t) -> ComplexDouble
        {
            if (vertices.empty())
            {
                return ComplexDouble(0.0, 0.0);
            }

            double normalizedT = std::fmod(t, 2.0 * M_PI);
            if (normalizedT < 0) normalizedT += 2.0 * M_PI;

            double indexF = normalizedT * vertices.size() / (2.0 * M_PI);
            int index1 = static_cast<int>(std::floor(indexF)) % vertices.size();
            int index2 = (index1 + 1) % vertices.size();
            double frac = indexF - index1;

            return vertices[index1] * ComplexDouble(1.0 - frac) + vertices[index2] * ComplexDouble(frac);
        };

        // Compute the derivative of the parameterization
        auto derivFunc = [vertices](double t) -> ComplexDouble
        {
            if (vertices.empty())
            {
                return ComplexDouble(0.0, 0.0);
            }

            double normalizedT = std::fmod(t, 2.0 * M_PI);
            if (normalizedT < 0) normalizedT += 2.0 * M_PI;

            double indexF = normalizedT * vertices.size() / (2.0 * M_PI);
            int index1 = static_cast<int>(std::floor(indexF)) % vertices.size();
            int index2 = (index1 + 1) % vertices.size();

            // Derivative is constant along each edge
            return (vertices[index2] - vertices[index1]) * ComplexDouble(vertices.size() / (2.0 * M_PI));
        };

        auto component = std::make_shared<AnalyticBoundaryComponent>(paramFunc, derivFunc);
        return std::make_shared<SimpleBoundary>(component);
    }
};

/**
 * @brief Multiply connected domain (connectivity >= 2)
 */
class MultiplyConnectedDomain : public Domain
{
protected:
    std::vector<std::shared_ptr<Boundary>> boundaries;

public:
    /**
     * @brief Construct a new Multiply Connected Domain
     * @param domainBoundaries Vector of boundaries
     * @param isExternalDomain Whether this is an external domain
     */
    MultiplyConnectedDomain(
        std::vector<std::shared_ptr<Boundary>> domainBoundaries,
        bool isExternalDomain = false)
        : Domain(isExternalDomain, domainBoundaries.size())
        , boundaries(domainBoundaries)
    {}

    void addBoundary(std::shared_ptr<Boundary> boundary)
    {
        boundaries.push_back(boundary);
        connectivity = boundaries.size();
    }

    std::vector<std::shared_ptr<Boundary>>& getBoundaries()
    {
        return boundaries;
    }

    const std::vector<std::shared_ptr<Boundary>>& getBoundaries() const
    {
        return boundaries;
    }

    bool contains(const ComplexDouble& z) const override
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
            auto samples = boundaries[i]->sample(500);
            if (samples.empty())
            {
                continue;
            }

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

    void transformBoundary(std::function<ComplexDouble(const ComplexDouble&)> transform) override
    {
        // Transform each boundary with increased sampling resolution (Key Change #3)
        const int transformSampleCount = 2000;

        for (auto& boundary : boundaries)
        {
            auto samples = boundary->sample(transformSampleCount);
            std::vector<ComplexDouble> transformedSamples;
            transformedSamples.reserve(samples.size());

            for (const auto& point : samples)
            {
                transformedSamples.push_back(transform(point));
            }

            // Create a new boundary from transformed samples
            auto discreteComponent = std::make_shared<DiscreteBoundaryComponent>(transformedSamples);
            boundary = std::make_shared<SimpleBoundary>(discreteComponent);
        }
    }
};

#endif // DOMAIN_HPP
