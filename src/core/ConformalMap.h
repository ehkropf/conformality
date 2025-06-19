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
#include "../domains/Domain.h"

#include <memory>

// Forward declaration to avoid circular dependency
class ConformalMapMethod;

/**
 * @brief Orchestrator class for conformal maps
 *
 * This class represents a conformal mapping between two domains and acts as
 * an orchestrator, coordinating between the domain geometry and the computational
 * method. It provides a clean interface for map evaluation while delegating
 * the actual computation to the method implementation.
 */
class ConformalMap
{
private:
    std::shared_ptr<Domain> mp_source_domain;
    std::shared_ptr<Domain> mp_target_domain;
    std::shared_ptr<ConformalMapMethod> mp_method;
    MappingType m_mapping_type;

public:
    /**
     * @brief Construct a new Conformal Map
     *
     * @param source Source domain
     * @param target Target domain
     * @param method_impl Method to use for computation
     */
    ConformalMap(
        std::shared_ptr<Domain> source,
        std::shared_ptr<Domain> target,
        std::shared_ptr<ConformalMapMethod> method_impl = nullptr
    );

    /**
     * @brief Copy constructor
     */
    ConformalMap(const ConformalMap& other) = default;

    /**
     * @brief Move constructor
     */
    ConformalMap(ConformalMap&& other) = default;

    /**
     * @brief Copy assignment operator
     */
    ConformalMap& operator=(const ConformalMap& other) = default;

    /**
     * @brief Move assignment operator
     */
    ConformalMap& operator=(ConformalMap&& other) = default;

    /**
     * @brief Default destructor
     */
    ~ConformalMap() = default;

    /**
     * @brief Map a point from the source domain to the target domain
     *
     * @param z Point in the source domain
     * @return Complex Mapped point in the target domain
     */
    Complex map(const Complex& z) const;

    /**
     * @brief Map a point from the target domain back to the source domain
     *
     * @param w Point in the target domain
     * @return Complex Mapped point in the source domain
     */
    Complex inverseMap(const Complex& w) const;

    /**
     * @brief Get the source domain
     *
     * @return std::shared_ptr<Domain> Source domain
     */
    std::shared_ptr<Domain> getSourceDomain() const
    {
        return mp_source_domain;
    }

    /**
     * @brief Get the target domain
     *
     * @return std::shared_ptr<Domain> Target domain
     */
    std::shared_ptr<Domain> getTargetDomain() const
    {
        return mp_target_domain;
    }

    /**
     * @brief Get the mapping type
     *
     * @return MappingType The type of this conformal mapping
     */
    MappingType getMappingType() const
    {
        return m_mapping_type;
    }

    /**
     * @brief Set the method used to compute the map
     *
     * @param method_impl Implementation of the method
     */
    void setMethod(std::shared_ptr<ConformalMapMethod> method_impl);

    /**
     * @brief Get the method used to compute the map
     *
     * @return std::shared_ptr<ConformalMapMethod> The method implementation
     */
    std::shared_ptr<ConformalMapMethod> getMethod() const
    {
        return mp_method;
    }

    /**
     * @brief Compute the map using the specified method
     *
     * @param target_accuracy Target accuracy for the computation
     */
    void compute(double target_accuracy = 1e-10);
};
