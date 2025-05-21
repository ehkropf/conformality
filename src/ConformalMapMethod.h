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

#ifndef CONFORMAL_MAP_METHOD_HPP
#define CONFORMAL_MAP_METHOD_HPP

#include <string>
#include "ConformalMap.h"

/**
 * @brief Abstract base class for conformal mapping methods
 *
 * This class defines the interface for methods that compute conformal mappings
 * between domains. Derived classes implement specific algorithms.
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
     * @brief Get the achieved accuracy of the last computation
     *
     * @return double Achieved accuracy
     */
    double getAccuracy() const;

    /**
     * @brief Get the number of iterations performed in the last computation
     *
     * @return int Number of iterations
     */
    int getIterationCount() const;

    /**
     * @brief Set the maximum number of iterations
     *
     * @param max Maximum number of iterations
     */
    void setMaxIterations(int max);

protected:
    /**
     * @brief Validate that the map is of the expected type
     *
     * @param map_instance Map to validate
     * @param expected_type Expected type name
     * @throws std::invalid_argument if the map is not of the expected type
     */
    void validateMapType(ConformalMap& map_instance, const std::string& expected_type) const;
};

#endif // CONFORMAL_MAP_METHOD_HPP
