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

#ifndef BOUNDARY_COMPONENT_HPP
#define BOUNDARY_COMPONENT_HPP

#include "Complex.h"

#include <functional>
#include <vector>

// Forward declarations
enum class InterpolationMethod {
    LINEAR,
    CUBIC_SPLINE,
    FOURIER
};

/**
 * @brief Abstract base class for boundary components
 */
class BoundaryComponent
{
protected:
    int index = 0;

public:
    virtual ~BoundaryComponent() = default;

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value in [0, 2π)
     * @return ComplexDouble point on the boundary
     */
    virtual ComplexDouble evaluate(double t) const = 0;

    /**
     * @brief Evaluate the derivative of the boundary at parameter t
     * @param t Parameter value in [0, 2π)
     * @return ComplexDouble derivative value
     */
    virtual ComplexDouble evaluateDerivative(double t) const = 0;

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    virtual std::vector<ComplexDouble> sample(size_t numPoints) const = 0;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z ComplexDouble point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    virtual double findParameterization(const ComplexDouble& z) const = 0;

    /**
     * @brief Get the index of this boundary component
     * @return Component index
     */
    int getIndex() const
    {
        return index;
    }

    /**
     * @brief Set the index of this boundary component
     * @param idx New index value
     */
    void setIndex(int idx)
    {
        index = idx;
    }
};

/**
 * @brief Boundary component defined by analytic functions
 */
class AnalyticBoundaryComponent : public BoundaryComponent
{
private:
    std::function<ComplexDouble(double)> parameterization;
    std::function<ComplexDouble(double)> derivative;

public:
    /**
     * @brief Construct a new Analytic Boundary Component
     * @param paramFunc Function mapping parameter t to complex point
     * @param derivFunc Function providing the derivative
     */
    AnalyticBoundaryComponent(std::function<ComplexDouble(double)> paramFunc,
                              std::function<ComplexDouble(double)> derivFunc);

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value
     * @return ComplexDouble point on the boundary
     */
    ComplexDouble evaluate(double t) const override
    {
        return parameterization(t);
    }

    /**
     * @brief Evaluate the derivative at parameter t
     * @param t Parameter value
     * @return ComplexDouble derivative
     */
    ComplexDouble evaluateDerivative(double t) const override
    {
        return derivative(t);
    }

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    std::vector<ComplexDouble> sample(size_t numPoints) const override;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z ComplexDouble point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    double findParameterization(const ComplexDouble& z) const override;
};

/**
 * @brief Boundary component defined by discrete points
 */
class DiscreteBoundaryComponent : public BoundaryComponent
{
private:
    std::vector<ComplexDouble> points;
    // FIXME: Uncomment this when implementing evaluation via method.
    // InterpolationMethod method;

public:
    /**
     * @brief Construct a new Discrete Boundary Component
     * @param pts Vector of points defining the boundary
     * @param interpolationMethod Method used for interpolation
     */
    DiscreteBoundaryComponent(const std::vector<ComplexDouble>& pts);

    // FIXME: Implementation method.
    //DiscreteBoundaryComponent(
    //    std::vector<ComplexDouble> pts,
    //    InterpolationMethod interpolationMethod = InterpolationMethod::LINEAR
    //)
    //    : points(pts)
    //    , method(interpolationMethod)
    //{}

    ComplexDouble evaluate(double t) const override;

    ComplexDouble evaluateDerivative(double t) const override;

    std::vector<ComplexDouble> sample(size_t numPoints) const override;

    double findParameterization(const ComplexDouble& z) const override;

    void resample(int numPoints);

    void smooth(double factor);
};

#endif // BOUNDARY_COMPONENT_HPP