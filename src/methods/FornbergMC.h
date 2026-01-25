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

#include "ConformalMapMethod.h"
#include "FornbergMCConfiguration.h"
#include "../core/Types.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

// Forward declarations
class MultiplyConnectedDomain;
class PMatrixBuilder;
class CGSolver;
class FornbergCanonicalDomain;

// Test helper forward declaration
#ifdef TESTING
#include <gtest/gtest_prod.h>
#endif

/**
 * @brief Implementation of Fornberg-like method for multiply connected domains
 *
 * This class implements the Fornberg-like method for numerical conformal mapping
 * of bounded multiply connected domains. The method uses Newton iteration to solve
 * for boundary correspondences S_ν(θ) and conformal moduli (c_ν, ρ_ν), employing
 * analyticity conditions on Fourier coefficients to enforce conformality.
 *
 * The algorithm computes conformal maps f: D → Ω where:
 * - Source domain D: Canonical domain (unit disk with circular holes)
 * - Target domain Ω: Bounded m-connected domain with smooth boundary components
 *
 * Key features:
 * - Newton iteration with adaptive damping
 * - Custom conjugate gradient solver for linear systems
 * - Automatic annulus detection and optimization (m=2)
 * - Comprehensive convergence monitoring
 * - Boundary point redistribution for numerical stability
 */
class FornbergMC : public ConformalMapMethod
{
#ifdef TESTING
    // Friend declaration for testing private methods
    FRIEND_TEST(FornbergMCFourierTest, CircularBoundaryCoefficients);
    FRIEND_TEST(FornbergMCFormSystemTest, DimensionsAnnulus);
    FRIEND_TEST(FornbergMCFormSystemTest, DimensionsGeneral);
    FRIEND_TEST(FornbergMCFormSystemTest, NonZeroOutput);
#endif

private:
    // Configuration
    FornbergMCConfiguration m_config;
    std::shared_ptr<MultiplyConnectedDomain> mp_user_domain{nullptr};        // Target domain (user's actual domain)
    std::shared_ptr<FornbergCanonicalDomain> mp_canonical_domain{nullptr};   // Source domain (unit disk + holes)

    // Matrix construction components
    std::unique_ptr<PMatrixBuilder> mp_matrix_builder{nullptr};
    std::unique_ptr<CGSolver> mp_cg_solver{nullptr};

    // Newton iteration state
    Eigen::MatrixXd m_S;                    // Boundary correspondences S_ν(θ)
    Eigen::VectorXcd m_conformal_moduli;    // Centers c_ν and radii ρ_ν

    // Linear system components
    Eigen::MatrixXcd m_D;                   // System matrix
    Eigen::VectorXcd m_g;                   // RHS vector
    Eigen::VectorXd m_U;                    // Solution vector

    // Series coefficients
    Eigen::MatrixXcd m_a;                   // Fourier coefficients a_{ν,j}

    // Algorithm state
    int m_connectivity{0};                     // Number of boundary components
    bool m_is_annulus{false};                  // Whether domain is an annulus (m=2)
    bool m_is_converged{false};                // Convergence flag
    double m_current_residual{std::numeric_limits<double>::max()};  // Current Newton update norm
    std::vector<double> m_residual_history; // Newton residual history

    // Auxiliary data for Newton updates
    Eigen::MatrixXd m_abs_eta;  // |eta| values (boundary derivative magnitudes), size: (m_connectivity, N)

    // Boundary parameterization
    std::vector<std::vector<Complex>> m_boundary_samples; // Sampled boundary points
    std::vector<std::vector<double>> m_parameter_values;  // Parameter values θ

public:
    /**
     * @brief Construct a new FornbergMC method
     * @param config Configuration parameters
     */
    explicit FornbergMC(const FornbergMCConfiguration& config = FornbergMCConfiguration{});

    /**
     * @brief Destructor
     */
    ~FornbergMC() override = default;

    /**
     * @brief Compute the conformal map using the Fornberg-like method
     * @param map_instance The map to compute
     * @param target_accuracy Target accuracy for convergence
     */
    void compute(ConformalMap& map_instance, double target_accuracy = 1e-6) override;

    /**
     * @brief Evaluate the computed map at a point
     * @param z Point in the canonical domain (source)
     * @return Complex Mapped point in the multiply connected domain (target)
     */
    Complex map(const Complex& z) const override;

    /**
     * @brief Evaluate the inverse of the computed map at a point
     * @param w Point in the multiply connected domain (target)
     * @return Complex Mapped point in the canonical domain (source)
     */
    Complex inverseMap(const Complex& w) const override;

    /**
     * @brief Get the current configuration
     * @return Configuration structure
     */
    const FornbergMCConfiguration& getConfiguration() const
    {
        return m_config;
    }

    /**
     * @brief Update the configuration (requires recomputation)
     * @param config New configuration
     */
    void setConfiguration(const FornbergMCConfiguration& config);

    /**
     * @brief Get the Fourier coefficients from the last computation
     * @return Matrix of coefficients a_{ν,j}
     */
    const Eigen::MatrixXcd& getFourierCoefficients() const
    {
        return m_a;
    }

    /**
     * @brief Get the conformal moduli from the last computation
     * @return Vector containing centers c_ν and radii ρ_ν
     */
    const Eigen::VectorXcd& getConformalModuli() const
    {
        return m_conformal_moduli;
    }

    /**
     * @brief Check if the last computation converged
     * @return True if converged
     */
    bool hasConverged() const
    {
        return m_is_converged;
    }

    /**
     * @brief Get the residual norm from the last computation
     * @return Current residual norm
     */
    double getCurrentResidual() const
    {
        return m_current_residual;
    }

    /**
     * @brief Get the full residual history
     * @return Vector of residual norms for each Newton iteration
     */
    const std::vector<double>& getResidualHistory() const
    {
        return m_residual_history;
    }

    /**
     * @brief Check if domain is detected as annulus case
     * @return True if annulus optimization is being used
     */
    bool isAnnulusCase() const
    {
        return m_is_annulus;
    }

    /**
     * @brief Get the system matrix D (for testing)
     * @return Reference to the D matrix
     */
    const Eigen::MatrixXcd& getSystemMatrix() const
    {
        return m_D;
    }

    /**
     * @brief Get the RHS vector g (for testing)
     * @return Reference to the g vector
     */
    const Eigen::VectorXcd& getRHSVector() const
    {
        return m_g;
    }

    /**
     * @brief Get the boundary correspondence matrix S (for testing)
     * @return Reference to the S matrix
     */
    const Eigen::MatrixXd& getBoundaryCorrespondences() const
    {
        return m_S;
    }

protected:
    /**
     * @brief Validate that the source domain is a canonical multiply connected domain
     * @param domain Source domain to validate
     * @throws std::invalid_argument if the domain is not compatible
     */
    void validateSourceDomain(std::shared_ptr<Domain> domain) const override;

    /**
     * @brief Validate that the target domain is a bounded multiply connected domain
     * @param domain Target domain to validate
     * @throws std::invalid_argument if the domain is not compatible
     */
    void validateTargetDomain(std::shared_ptr<Domain> domain) const override;

private:
    /**
     * @brief Initialize the Newton iteration framework
     */
    void initializeNewtonIteration();

    /**
     * @brief Form the linear system D*U = g
     * Detects annulus case and uses appropriate formulation
     */
    void formSystem();

    /**
     * @brief Solve the linear system using conjugate gradient method
     */
    void solveSystem();

    /**
     * @brief Perform Newton update of boundary correspondences and conformal moduli
     */
    void newtonUpdate();

    /**
     * @brief Check convergence based on Newton update norm
     * @param tolerance Convergence tolerance to use
     * @return True if converged
     */
    bool checkConvergence(double tolerance);

    /**
     * @brief Detect if domain configuration represents an annulus (m=2)
     * @return True if annulus case detected
     */
    bool detectAnnulusCase() const;

    /**
     * @brief Initialize conformal moduli with geometric heuristics
     */
    void initializeConformalModuli();

    /**
     * @brief Sample boundary components at parameter values
     */
    void sampleBoundaries();

    /**
     * @brief Redistribute boundary parameters if needed
     * @return True if redistribution was performed
     */
    bool redistributeBoundaryParameters();

    /**
     * @brief Compute Fourier coefficients from current boundary data
     */
    void computeFourierCoefficients();

    /**
     * @brief Validate that N is a power of 2 for FFT compatibility
     * @param N Number of points to validate
     * @return True if valid
     */
    bool isPowerOfTwo(int N) const;

    /**
     * @brief Print iteration diagnostics if verbose mode enabled
     * @param iteration Current iteration number
     */
    void printIterationDiagnostics(size_t iteration) const;
};
