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

#include "Types.h"
#include "StatusManager.h"

#include <functional>
#include <vector>
#include <memory>

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
    std::shared_ptr<IStatusManager> p_statusManager;

public:
    virtual ~BoundaryComponent() = default;

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value in [0, 2π)
     * @return Complex point on the boundary
     */
    virtual Complex evaluate(double t) const = 0;

    /**
     * @brief Evaluate the derivative of the boundary at parameter t
     * @param t Parameter value in [0, 2π)
     * @return Complex derivative value
     */
    virtual Complex evaluateDerivative(double t) const = 0;

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    virtual std::vector<Complex> sample(size_t numPoints) const = 0;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z Complex point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    virtual double findParameterization(const Complex& z) const = 0;


    /**
     * @brief Set the status manager for this boundary component
     * @param manager Shared pointer to status manager
     */
    void setStatusManager(std::shared_ptr<IStatusManager> manager)
    {
        p_statusManager = manager;
    }
};

/**
 * @brief Boundary component defined by analytic functions
 */
class AnalyticBoundaryComponent : public BoundaryComponent
{
private:
    std::function<Complex(double)> parameterization;
    std::function<Complex(double)> derivative;

public:
    /**
     * @brief Construct a new Analytic Boundary Component
     * @param paramFunc Function mapping parameter t to complex point
     * @param derivFunc Function providing the derivative
     * @param statusMgr Optional status manager for reporting warnings/errors
     */
    AnalyticBoundaryComponent(std::function<Complex(double)> paramFunc,
                              std::function<Complex(double)> derivFunc,
                              std::shared_ptr<IStatusManager> statusMgr = nullptr);

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value
     * @return Complex point on the boundary
     */
    Complex evaluate(double t) const override
    {
        return parameterization(t);
    }

    /**
     * @brief Evaluate the derivative at parameter t
     * @param t Parameter value
     * @return Complex derivative
     */
    Complex evaluateDerivative(double t) const override
    {
        return derivative(t);
    }

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    std::vector<Complex> sample(size_t numPoints) const override;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z Complex point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    double findParameterization(const Complex& z) const override;
};

/**
 * @brief Boundary component defined by discrete points
 */
class DiscreteBoundaryComponent : public BoundaryComponent
{
private:
    std::vector<Complex> points;
    InterpolationMethod method;

public:
    /**
     * @brief Construct a new Discrete Boundary Component
     * @param pts Vector of points defining the boundary
     * @param interpolationMethod Method used for interpolation
     * @param statusMgr Optional status manager for reporting warnings/errors
     */
    DiscreteBoundaryComponent(
            const std::vector<Complex>& pts,
            InterpolationMethod interpolationMethod = InterpolationMethod::LINEAR,
            std::shared_ptr<IStatusManager> statusMgr = nullptr);

    Complex evaluate(double t) const override;

    Complex evaluateDerivative(double t) const override;

    std::vector<Complex> sample(size_t numPoints) const override;

    double findParameterization(const Complex& z) const override;

    void resample(int numPoints);

    void smooth(double factor);
};

#endif // BOUNDARY_COMPONENT_HPP
