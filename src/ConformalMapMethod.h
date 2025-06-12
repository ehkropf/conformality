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

#include "Types.h"

#include <memory>

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
    double m_achieved_accuracy;
    int m_max_iterations;
    int m_iteration_count;

public:
    /**
     * @brief Construct a new Conformal Map Method with default settings
     */
    ConformalMapMethod();

    /**
     * @brief Construct a new Conformal Map Method with custom settings
     * @param max_iter Maximum number of iterations (must be positive)
     */
    explicit ConformalMapMethod(int max_iter);

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
        return m_achieved_accuracy;
    }

    /**
     * @brief Get the number of iterations performed in the last computation
     *
     * @return int Number of iterations
     */
    int getIterationCount() const
    {
        return m_iteration_count;
    }

    /**
     * @brief Set the maximum number of iterations
     *
     * @param max Maximum number of iterations
     */
    void setMaxIterations(int max);

    /**
     * @brief Validate a domain for use with this method
     *
     * Performs both compatibility and geometry validation.
     *
     * @param domain Domain to validate
     * @param expected_connectivity Expected connectivity (0 = simply connected, etc.)
     * @throws std::invalid_argument if the domain is not valid for this method
     */
    void validateDomain(std::shared_ptr<Domain> domain, int expected_connectivity) const;

    /**
     * @brief Validate both source and target domains for use with this method
     *
     * Validates that both domains in the conformal map are compatible with this method.
     * Derived classes can override this to specify domain-specific connectivity requirements.
     *
     * @param map_instance The conformal map whose domains should be validated
     * @throws std::invalid_argument if either domain is not valid for this method
     */
    virtual void validateDomains(const ConformalMap& map_instance) const;

protected:
    /**
     * @brief Validate that the domain is compatible with this method
     *
     * Checks connectivity and null pointer. Always called by validateDomain.
     *
     * @param domain Domain to validate
     * @param expected_connectivity Expected connectivity (0 = simply connected, etc.)
     * @throws std::invalid_argument if the domain is not compatible
     */
    void validateDomainCompatibility(std::shared_ptr<Domain> domain, int expected_connectivity) const;

    /**
     * @brief Validate that the domain has the required geometric properties
     *
     * Method-specific geometric validation. Must be implemented by derived classes.
     *
     * @param domain Domain to validate
     * @throws std::invalid_argument if the domain doesn't have required properties
     */
    virtual void validateDomainGeometry(std::shared_ptr<Domain> domain) const = 0;
};

