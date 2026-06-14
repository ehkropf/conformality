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

#include "../core/MethodInfo.h"
#include "../core/Types.h"

#include <functional>
#include <memory>

// Forward declarations to avoid circular dependency
class ConformalMap;
class Domain;

/**
 * @brief Typed progress data from iterative methods
 *
 * Provides machine-readable progress information separate from human-readable
 * StatusMessage logging. GUI components consume this directly without parsing strings.
 */
struct ProgressUpdate
{
    int iteration{0};       ///< 1-based iteration number
    double residual{0.0};   ///< Current residual/error metric
};

/**
 * @brief Abstract base class for conformal mapping methods
 *
 * This class defines the interface for methods that compute conformal mappings
 * between domains. Methods own all algorithm-specific state and handle both
 * computation and map evaluation. They validate domain compatibility at runtime.
 */
class ConformalMapMethod
{
public:
    enum class DomainRole
    {
        Source,
        Target
    };

    using ProgressCallback = std::function<void(const ProgressUpdate&)>;

protected:
    double m_achieved_accuracy;
    int m_max_iterations;
    int m_iteration_count;
    std::function<bool()> m_cancellationCheck;  ///< Returns true if computation should be cancelled
    ProgressCallback m_progressCallback;        ///< Optional typed progress reporting callback

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
     * @brief Get a descriptor with method name, parameters, and results
     * @return conformality::MethodInfo populated by the concrete method (empty by default)
     */
    virtual conformality::MethodInfo getMethodInfo() const { return {}; }

    /**
     * @brief Set a cancellation check callback for cooperative cancellation
     * @param check Function returning true if computation should be cancelled, or nullptr to clear
     */
    void setCancellationCheck(std::function<bool()> check)
    {
        m_cancellationCheck = std::move(check);
    }

    /**
     * @brief Set a typed progress callback for live iteration reporting
     *
     * Unlike StatusMessage logging (human-readable text), this provides typed
     * fields for GUI consumption without string parsing.
     *
     * @param callback Function to invoke with each progress update, or nullptr to clear
     */
    void setProgressCallback(ProgressCallback callback)
    {
        m_progressCallback = std::move(callback);
    }

    /**
     * @brief Validate a domain for use with this method
     *
     * Performs both compatibility and geometry validation.
     *
     * @param domain Domain to validate
     * @param expected_connectivity Expected connectivity (0 = simply connected, etc.)
     * @param role Whether validating source or target domain
     * @throws std::invalid_argument if the domain is not valid for this method
     */
    void validateDomain(std::shared_ptr<Domain> domain, int expected_connectivity, DomainRole role) const;

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
     * @brief Validate that the source domain has the required geometric properties
     *
     * Method-specific validation for source domain geometry. Must be implemented by derived classes.
     *
     * @param domain Source domain to validate
     * @throws std::invalid_argument if the domain doesn't have required properties for source
     */
    virtual void validateSourceDomain(std::shared_ptr<Domain> domain) const = 0;

    /**
     * @brief Validate that the target domain has the required geometric properties
     *
     * Method-specific validation for target domain geometry. Must be implemented by derived classes.
     *
     * @param domain Target domain to validate
     * @throws std::invalid_argument if the domain doesn't have required properties for target
     */
    virtual void validateTargetDomain(std::shared_ptr<Domain> domain) const = 0;
};
