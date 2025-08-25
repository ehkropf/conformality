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

#include <vector>
#include <concepts>
#include <cassert>

/**
 * @brief Concept for types that can be used with polyval
 * 
 * Types must support multiplication, addition, and copy construction.
 * This includes arithmetic types and complex number types.
 */
template<typename T>
concept PolyvalCompatible = std::is_arithmetic_v<T> || requires(T a, T b) {
    { a * b } -> std::convertible_to<T>;
    { a + b } -> std::convertible_to<T>;
    T{a}; // Copy constructible
};

/**
 * @brief Evaluate polynomial with coefficients in descending power order
 * 
 * This function evaluates a polynomial using Horner's method for efficiency.
 * The coefficients are expected in descending power order: [c_n, c_{n-1}, ..., c_1, c_0]
 * to represent the polynomial c_n*x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0.
 * 
 * @tparam T Coefficient and evaluation point type (must satisfy PolyvalCompatible)
 * @param coeffs Vector of coefficients in descending power order (must not be empty)
 * @param x Point at which to evaluate the polynomial
 * @return T Polynomial value at x
 * 
 * @pre coeffs is not empty
 */
template<PolyvalCompatible T>
T polyval(const std::vector<T>& coeffs, const T& x)
{
    // Contract: precondition check
    assert(!coeffs.empty() && "polyval: coefficient vector cannot be empty");
    
    // Horner's method: evaluate polynomial efficiently in O(n) operations
    // For coeffs = [c_n, c_{n-1}, ..., c_1, c_0], compute c_n*x^n + ... + c_0
    // Algorithm: result = c_n, then result = result*x + c_{n-1}, etc.
    
    T result = coeffs[0]; // Start with highest degree coefficient
    
    for (size_t i = 1; i < coeffs.size(); ++i)
    {
        result = result * x + coeffs[i]; // result*x + next_coeff
    }
    
    return result;
}

