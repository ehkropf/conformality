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

#pragma once

#include "BoundaryComponent.h"

#include <memory>
#include <vector>

/**
 * @brief Boundary consisting of one or more boundary components
 *
 * Handles both simply-connected (single component) and multiply-connected
 * (multiple components) boundaries automatically.
 */
class Boundary
{
protected:
    std::vector<std::shared_ptr<BoundaryComponent>> components;

public:
    Boundary() = default;

    /**
     * @brief Construct boundary with components
     * @param boundaryComponents Vector of boundary components
     */
    Boundary(std::vector<std::shared_ptr<BoundaryComponent>> boundaryComponents);

    /**
     * @brief Construct boundary with single component
     * @param component Single boundary component
     */
    Boundary(std::shared_ptr<BoundaryComponent> component);

    /**
     * @brief Evaluate the boundary at parameter t on specified component
     * @param t Parameter value in [0, totalLength())
     * @param componentIndex Index of the boundary component (default: 0)
     * @return Complex point on the boundary
     */
    Complex evaluate(double t, size_t componentIndex = 0) const;

    /**
     * @brief Evaluate the derivative at parameter t on specified component
     * @param t Parameter value in [0, totalLength())
     * @param componentIndex Index of the boundary component (default: 0)
     * @return Complex derivative
     */
    Complex evaluateDerivative(double t, size_t componentIndex = 0) const;

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample from each component
     * @return Vector of vectors, one for each boundary component
     */
    std::vector<std::vector<Complex>> sample(size_t numPoints) const;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z Complex point on or near the boundary
     * @param componentIndex Index of the boundary component to search on
     * @return Parameter value t such that evaluate(t, componentIndex) ≈ z
     */
    double findParameterization(const Complex& z, size_t componentIndex = 0) const;

    /**
     * @brief Get the total parameter length for a boundary component
     * @param componentIndex Index of the boundary component (default: 0)
     * @return Total parameter length (delegates to component's totalLength())
     */
    double totalLength(size_t componentIndex = 0) const;

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
