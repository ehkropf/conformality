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

#ifndef CONFORMAL_MAP_HPP
#define CONFORMAL_MAP_HPP

#include <memory>
#include "Types.h"
#include "Domain.h"

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
    std::shared_ptr<Domain> source_domain;
    std::shared_ptr<Domain> target_domain;
    bool is_external;
    std::shared_ptr<ConformalMapMethod> method;

public:
    /**
     * @brief Construct a new Conformal Map
     *
     * @param source Source domain
     * @param target Target domain
     * @param method_impl Method to use for computation
     * @param external Whether this is an external map
     */
    ConformalMap(
        std::shared_ptr<Domain> source,
        std::shared_ptr<Domain> target,
        std::shared_ptr<ConformalMapMethod> method_impl = nullptr,
        bool external = false
    );

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
        return source_domain;
    }

    /**
     * @brief Get the target domain
     *
     * @return std::shared_ptr<Domain> Target domain
     */
    std::shared_ptr<Domain> getTargetDomain() const
    {
        return target_domain;
    }

    /**
     * @brief Check if this is an external map
     *
     * @return true if external, false if internal
     */
    bool isExternalMap() const
    {
        return is_external;
    }

    /**
     * @brief Set whether this is an external map
     *
     * @param external True for external, false for internal
     */
    void setExternal(bool external)
    {
        is_external = external;
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
        return method;
    }

    /**
     * @brief Compute the map using the specified method
     *
     * @param target_accuracy Target accuracy for the computation
     */
    void compute(double target_accuracy = 1e-10);
};

#endif // CONFORMAL_MAP_HPP
