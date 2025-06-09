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

#ifndef CONFORMAL_MAP_METHOD_HPP
#define CONFORMAL_MAP_METHOD_HPP

#include "Types.h"

#include <memory>
#include <string>

// Forward declarations to avoid circular dependency
class ConformalMap;
class Domain;

/**
 * @brief Abstract base class for conformal mapping methods
 *
 * This class defines the interface for methods that compute conformal mappings
 * between domains. Methods own all algorithm-specific state and handle both
 * computation and map evaluation. They validate domain compatibility at runtime.
 */
class ConformalMapMethod
{
protected:
    double accuracy;
    int max_iterations;
    int iteration_count;

public:
    /**
     * @brief Construct a new Conformal Map Method
     */
    ConformalMapMethod();

    /**
     * @brief Virtual destructor
     */
    virtual ~ConformalMapMethod() = default;

    /**
     * @brief Compute the conformal map
     *
     * @param map_instance The map to compute
     * @param target_accuracy Target accuracy for the computation
     */
    virtual void compute(ConformalMap& map_instance, double target_accuracy = 1e-10) = 0;

    /**
     * @brief Evaluate the computed map at a point
     *
     * @param z Point in the source domain
     * @return Complex Mapped point in the target domain
     */
    virtual Complex map(const Complex& z) const = 0;

    /**
     * @brief Evaluate the inverse of the computed map at a point
     *
     * @param w Point in the target domain
     * @return Complex Mapped point in the source domain
     */
    virtual Complex inverseMap(const Complex& w) const = 0;

    /**
     * @brief Get the achieved accuracy of the last computation
     *
     * @return double Achieved accuracy
     */
    double getAccuracy() const
    {
        return accuracy;
    }

    /**
     * @brief Get the number of iterations performed in the last computation
     *
     * @return int Number of iterations
     */
    int getIterationCount() const
    {
        return iteration_count;
    }

    /**
     * @brief Set the maximum number of iterations
     *
     * @param max Maximum number of iterations
     */
    void setMaxIterations(int max);

protected:
    /**
     * @brief Validate that the domain is compatible with this method
     *
     * @param domain Domain to validate
     * @param expected_connectivity Expected connectivity (0 = simply connected, etc.)
     * @throws std::invalid_argument if the domain is not compatible
     */
    void validateDomainCompatibility(std::shared_ptr<Domain> domain, int expected_connectivity) const;

    /**
     * @brief Validate that the domain has the required geometric properties
     *
     * @param domain Domain to validate
     * @throws std::invalid_argument if the domain doesn't have required properties
     */
    void validateDomainGeometry(std::shared_ptr<Domain> domain) const;
};

#endif // CONFORMAL_MAP_METHOD_HPP
