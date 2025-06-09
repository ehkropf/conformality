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

#include <complex>
#include <cmath>

/**
 * @brief A wrapper class for complex numbers.
 *
 * This class provides a consistent interface for complex numbers,
 * defaulting to std::complex<double> but supporting other complex
 * types through templating.
 *
 * @tparam T The underlying complex type (default: std::complex<double>)
 */
template <typename T = std::complex<double>>
class Complex
{
private:
    T value;

public:
    /**
     * @brief Default constructor.
     */
    Complex() : value(0.0, 0.0) {}

    /**
     * @brief Construct from real and imaginary parts.
     *
     * @param real The real part
     * @param imag The imaginary part
     */
    Complex(double real, double imag = 0.0) : value(real, imag) {}

    /**
     * @brief Construct from underlying complex type.
     *
     * @param complex The complex value to wrap
     */
    explicit Complex(const T& complex) : value(complex) {}

    /**
     * @brief Get the real part.
     *
     * @return The real part
     */
    double real() const
    {
        return std::real(value);
    }

    /**
     * @brief Get the imaginary part.
     *
     * @return The imaginary part
     */
    double imag() const
    {
        return std::imag(value);
    }

    /**
     * @brief Get the modulus (absolute value).
     *
     * @return The modulus
     */
    double abs() const
    {
        return std::abs(value);
    }

    /**
     * @brief Get the argument (phase angle).
     *
     * @return The argument in radians
     */
    double arg() const
    {
        return std::arg(value);
    }

    /**
     * @brief Get the underlying complex value.
     *
     * @return The underlying complex value
     */
// FIXME:  
    const T& getValue const
    {
        return value;
    }

    /**
     * @brief Get the underlying complex value.
     *
     * @return The underlying complex value
     */
    T& getValue()
    {
        return value;
    }

    /**
     * @brief Addition operator.
     *
     * @param other The complex number to add
     * @return The sum
     */
    Complex operator+(const Complex& other) const
    {
        return Complex(value + other.value);
    }

    /**
     * @brief Subtraction operator.
     *
     * @param other The complex number to subtract
     * @return The difference
     */
    Complex operator-(const Complex& other) const
    {
        return Complex(value - other.value);
    }

    /**
     * @brief Multiplication operator.
     *
     * @param other The complex number to multiply by
     * @return The product
     */
    Complex operator*(const Complex& other) const
    {
        return Complex(value * other.value);
    }

    /**
     * @brief Division operator.
     *
     * @param other The complex number to divide by
     * @return The quotient
     */
    Complex operator/(const Complex& other) const
    {
        return Complex(value / other.value);
    }

    /**
     * @brief Equality operator.
     *
     * @param other The complex number to compare with
     * @return true if the two complex numbers are equal, false otherwise
     */
    bool operator==(const Complex& other) const
    {
        return value == other.value;
    }

    /**
     * @brief Create a complex number from polar coordinates.
     *
     * @param modulus The modulus (absolute value)
     * @param argument The argument (phase angle) in radians
     * @return The complex number
     */
    static Complex fromPolar(double modulus, double argument)
    {
        return Complex(std::polar(modulus, argument));
    }
};

// Type alias for the default Complex type
using ComplexDouble = Complex<std::complex<double>>;
