/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
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
#include "Complex.hpp"
#include "Domain.hpp"

// Forward declaration to avoid circular dependency
class ConformalMapMethod;

/**
 * @brief Abstract base class for conformal maps
 *
 * This class represents a conformal mapping between two domains. It provides
 * the interface for evaluating the map and its inverse, as well as access to
 * the source and target domains.
 */
class ConformalMap
{
protected:
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
     * @param external Whether this is an external map
     */
    ConformalMap(
        std::shared_ptr<Domain> source,
        std::shared_ptr<Domain> target,
        bool external = false
    );

    /**
     * @brief Virtual destructor
     */
    virtual ~ConformalMap() = default;

    /**
     * @brief Map a point from the source domain to the target domain
     *
     * @param z Point in the source domain
     * @return ComplexDouble Mapped point in the target domain
     */
    virtual ComplexDouble map(const ComplexDouble& z) const = 0;

    /**
     * @brief Map a point from the target domain back to the source domain
     *
     * @param w Point in the target domain
     * @return ComplexDouble Mapped point in the source domain
     */
    virtual ComplexDouble inverseMap(const ComplexDouble& w) const = 0;

    /**
     * @brief Get the source domain
     *
     * @return std::shared_ptr<Domain> Source domain
     */
    std::shared_ptr<Domain> getSourceDomain() const;

    /**
     * @brief Get the target domain
     *
     * @return std::shared_ptr<Domain> Target domain
     */
    std::shared_ptr<Domain> getTargetDomain() const;

    /**
     * @brief Check if this is an external map
     *
     * @return true if external, false if internal
     */
    bool isExternalMap() const;

    /**
     * @brief Set whether this is an external map
     *
     * @param external True for external, false for internal
     */
    void setExternal(bool external);

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
    std::shared_ptr<ConformalMapMethod> getMethod() const;

    /**
     * @brief Compute the map using the specified method
     *
     * @param target_accuracy Target accuracy for the computation
     */
    void compute(double target_accuracy = 1e-10);
};

#endif // CONFORMAL_MAP_HPP
