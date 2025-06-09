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

#ifndef THEODORSEN_METHOD_HPP
#define THEODORSEN_METHOD_HPP

#include "ConformalMapMethod.h"
#include "FFTWWrapper.h"
#include <vector>

/**
 * @brief Implementation of Theodorsen's method for conformal mapping
 *
 * Theodorsen's method computes conformal maps from starlike domains to the unit disk
 * using an iterative FFT-based algorithm. This implementation supports both internal
 * and external mappings.
 *
 * The method works by iteratively refining a boundary correspondence function using
 * the conjugation operator and FFTs until convergence is achieved.
 */
class TheodorsenMethod : public ConformalMapMethod
{
private:
    size_t n_points;                              // Number of sample points (must be power of 2)
    std::vector<Complex> laurent_coeffs;    // Laurent series coefficients
    std::vector<Complex> boundary_samples;  // Sampled boundary points
    std::vector<double> phi_sequence;             // Angle sequence for sampling
    bool is_converged;                            // Convergence flag
    double residual_norm;                         // Last residual norm

    /**
     * @brief Sample the boundary of a starlike domain
     * @param domain Starlike domain to sample
     * @param external Whether this is an external mapping
     * @return Vector of boundary sample points
     */
    std::vector<Complex> sampleBoundary(std::shared_ptr<Domain> domain, bool external) const;

    /**
     * @brief Compute the boundary correspondence function rho(phi)
     * @param domain Starlike domain
     * @param external Whether this is an external mapping
     * @return Vector of boundary moduli
     */
    std::vector<double> computeBoundaryModuli(std::shared_ptr<Domain> domain, bool external) const;

    /**
     * @brief Perform one iteration of Theodorsen's algorithm
     * @param rho Current boundary moduli
     * @param external Whether this is an external mapping
     * @return Updated boundary moduli
     */
    std::vector<double> theodorsenIteration(const std::vector<double>& rho, bool external);

    /**
     * @brief Compute Laurent series coefficients from boundary data
     * @param boundary_data Boundary correspondence function values
     */
    void computeLaurentCoefficients(const std::vector<Complex>& boundary_data);

    /**
     * @brief Check convergence based on residual norm
     * @param old_rho Previous boundary moduli
     * @param new_rho Current boundary moduli
     * @param tolerance Convergence tolerance
     * @return True if converged
     */
    bool checkConvergence(const std::vector<double>& old_rho, const std::vector<double>& new_rho, double tolerance);

public:
    /**
     * @brief Construct a new Theodorsen Method
     * @param num_points Number of sample points (must be power of 2)
     */
    explicit TheodorsenMethod(size_t num_points = 128);

    /**
     * @brief Destructor
     */
    ~TheodorsenMethod() override = default;

    /**
     * @brief Compute the conformal map using Theodorsen's method
     * @param map_instance The map to compute
     * @param target_accuracy Target accuracy for convergence
     */
    void compute(ConformalMap& map_instance, double target_accuracy = 1e-10) override;

    /**
     * @brief Evaluate the computed map at a point
     * @param z Point in the source domain
     * @return Complex Mapped point in the target domain
     */
    Complex map(const Complex& z) const override;

    /**
     * @brief Evaluate the inverse of the computed map at a point
     * @param w Point in the target domain
     * @return Complex Mapped point in the source domain
     */
    Complex inverseMap(const Complex& w) const override;

    /**
     * @brief Get the number of sample points
     * @return Number of sample points
     */
    size_t getNumPoints() const
    {
        return n_points;
    }

    /**
     * @brief Set the number of sample points
     * @param num_points New number of sample points (must be power of 2)
     */
    void setNumPoints(size_t num_points);

    /**
     * @brief Get the Laurent series coefficients
     * @return Vector of Laurent coefficients
     */
    const std::vector<Complex>& getLaurentCoefficients() const
    {
        return laurent_coeffs;
    }

    /**
     * @brief Check if the last computation converged
     * @return True if converged
     */
    bool hasConverged() const
    {
        return is_converged;
    }

    /**
     * @brief Get the residual norm from the last computation
     * @return Residual norm
     */
    double getResidualNorm() const
    {
        return residual_norm;
    }

    /**
     * @brief Get the boundary sample points from the last computation
     * @return Vector of boundary sample points
     */
    const std::vector<Complex>& getBoundarySamples() const
    {
        return boundary_samples;
    }

private:
    /**
     * @brief Validate that the number of points is a power of 2
     * @param num_points Number to validate
     * @return True if valid
     */
    bool isPowerOfTwo(size_t num_points) const;
};

#endif // THEODORSEN_METHOD_HPP