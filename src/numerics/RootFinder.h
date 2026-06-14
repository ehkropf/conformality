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

#ifndef ROOT_FINDER_HPP
#define ROOT_FINDER_HPP

#include <functional>
#include <stdexcept>
#include <complex>

#include "../core/Tolerances.h"

/**
 * @brief Utility class for root finding and optimization algorithms
 */
class RootFinder
{
public:
    /**
     * @brief Find minimum of a scalar function using ternary search
     * @param objective Function to minimize (takes double, returns double)
     * @param low Lower bound of search interval
     * @param high Upper bound of search interval
     * @param tolerance Convergence tolerance
     * @param maxIterations Maximum number of iterations
     * @return Parameter value that minimizes the objective function
     */
    static double ternarySearch(std::function<double(double)> objective,
                                double low, double high,
                                double tolerance = 1e-9,
                                int maxIterations = 100);

    /**
     * @brief Newton's method for finding roots
     * @tparam T Function domain/range type (double or std::complex<double>)
     * @param function Function whose root to find
     * @param derivative Derivative of the function
     * @param initial Initial guess
     * @param tolerance Convergence tolerance
     * @param maxIterations Maximum number of iterations
     * @return Root of the function
     */
    template<typename T>
    static T newton(std::function<T(T)> function,
                    std::function<T(T)> derivative,
                    T initial,
                    double tolerance = 1e-9,
                    int maxIterations = 100);

    /**
     * @brief Exception thrown when root finding fails to converge
     */
    class ConvergenceError : public std::runtime_error
    {
    public:
        explicit ConvergenceError(const std::string& message)
            : std::runtime_error(message)
        {
        }
    };
};

// Template implementations
#include <cmath>
#include <algorithm>

template<typename T>
T RootFinder::newton(std::function<T(T)> function,
                     std::function<T(T)> derivative,
                     T initial,
                     double tolerance,
                     int maxIterations)
{
    T x = initial;

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        T fx = function(x);
        T dfx = derivative(x);

        if (std::abs(fx) < tolerance)
        {
            return x;
        }

        if (std::abs(dfx) < PIVOT_EPS)
        {
            throw ConvergenceError("Newton's method: derivative too small");
        }

        T newX = x - fx / dfx;

        double stepSize = std::abs(newX - x);
        double scale = std::max(std::abs(newX), std::abs(x));
        if (stepSize < tolerance * std::max(scale, 1.0))
        {
            return newX;
        }

        x = newX;
    }

    throw ConvergenceError("Newton's method failed to converge within " +
                           std::to_string(maxIterations) + " iterations");
}

#endif // ROOT_FINDER_HPP
