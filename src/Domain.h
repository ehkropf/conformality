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

#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Boundary.h"
#include <functional>
#include <memory>
#include <vector>

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
    int calculateWindingNumber(const ComplexDouble& z, const std::vector<ComplexDouble>& samples, double tolerance = 1e-12) const;

    /**
     * @brief Calculate distance from point to line segment
     * @param point Point to measure distance from
     * @param segmentStart Start of line segment
     * @param segmentEnd End of line segment
     * @return Distance to closest point on segment
     */
    double distanceToLineSegment(const ComplexDouble& point, const ComplexDouble& segmentStart, const ComplexDouble& segmentEnd) const;

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

    bool contains(const ComplexDouble& z) const override;

    void transformBoundary(std::function<ComplexDouble(const ComplexDouble&)> transform) override;
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

    bool contains(const ComplexDouble& z) const override;

private:
    static std::shared_ptr<Boundary> createBoundary(
        const ComplexDouble& center,
        std::function<double(double)> radiusFunc
    );
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

    bool contains(const ComplexDouble& z) const override;
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
    void setVertices(const std::vector<ComplexDouble>& newVertices);

private:
    static std::shared_ptr<Boundary> createBoundary(
        const std::vector<ComplexDouble>& vertices
    );
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

    void addBoundary(std::shared_ptr<Boundary> boundary);

    std::vector<std::shared_ptr<Boundary>>& getBoundaries()
    {
        return boundaries;
    }
    
    const std::vector<std::shared_ptr<Boundary>>& getBoundaries() const
    {
        return boundaries;
    }

    bool contains(const ComplexDouble& z) const override;

    void transformBoundary(std::function<ComplexDouble(const ComplexDouble&)> transform) override;
};

#endif // DOMAIN_HPP
