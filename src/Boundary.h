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

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "BoundaryComponent.h"

#include <memory>
#include <vector>

/**
 * @brief Abstract base class for boundaries
 *
 * A boundary consists of one or more boundary components
 */
class Boundary
{
protected:
    std::vector<std::shared_ptr<BoundaryComponent>> components;

public:
    virtual ~Boundary() = default;

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value
     * @return ComplexDouble point on the boundary
     */
    virtual ComplexDouble evaluate(double t) const = 0;

    /**
     * @brief Evaluate the derivative at parameter t
     * @param t Parameter value
     * @return ComplexDouble derivative
     */
    virtual ComplexDouble evaluateDerivative(double t) const = 0;

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    virtual std::vector<ComplexDouble> sample(int numPoints) const = 0;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z ComplexDouble point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    virtual double findParameterization(const ComplexDouble& z) const = 0;

    /**
     * @brief Add a boundary component
     * @param component Boundary component to add
     */
    void addComponent(std::shared_ptr<BoundaryComponent> component);

    /**
     * @brief Get the number of boundary components
     * @return Number of components
     */
    size_t getNumComponents() const
    {
        return components.size();
    }

    /**
     * @brief Get a specific boundary component
     * @param index Component index
     * @return Boundary component
     */
    const BoundaryComponent& getComponent(size_t index) const;

    /**
     * @brief Get all boundary components
     * @return Vector of boundary components
     */
    const std::vector<std::shared_ptr<BoundaryComponent>>& getComponents() const
    {
        return components;
    }
};

/**
 * @brief Simple boundary consisting of a single component
 */
class SimpleBoundary : public Boundary
{
public:
    /**
     * @brief Construct a new Simple Boundary
     * @param component Single boundary component
     */
    SimpleBoundary(std::shared_ptr<BoundaryComponent> component);

    ComplexDouble evaluate(double t) const override
    {
        return components[0]->evaluate(t);
    }

    ComplexDouble evaluateDerivative(double t) const override
    {
        return components[0]->evaluateDerivative(t);
    }

    std::vector<ComplexDouble> sample(int numPoints) const override
    {
        return components[0]->sample(numPoints);
    }

    double findParameterization(const ComplexDouble& z) const override
    {
        return components[0]->findParameterization(z);
    }
};

/**
 * @brief Composite boundary consisting of multiple components
 */
class CompositeBoundary : public Boundary
{
private:
    std::vector<double> parameterRanges;  // End of parameter range for each component

public:
    CompositeBoundary() = default;

    /**
     * @brief Construct a new Composite Boundary
     * @param components Vector of boundary components
     */
    CompositeBoundary(std::vector<std::shared_ptr<BoundaryComponent>> boundaryComponents);

    void addComponent(std::shared_ptr<BoundaryComponent> component);

    ComplexDouble evaluate(double t) const override;

    ComplexDouble evaluateDerivative(double t) const override;

    std::vector<ComplexDouble> sample(int numPoints) const override;

    double findParameterization(const ComplexDouble& z) const override;

private:
    void updateParameterRanges();
};

#endif // BOUNDARY_HPP