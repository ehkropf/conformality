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

#include "../core/Types.h"
#include "Domain.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

/**
 * @brief Canonical domain for the Fornberg method (unit disk with circular holes)
 *
 * This class represents the canonical source domain used in the Fornberg-like method:
 * a unit disk with (m-1) circular holes for m-connected domains, or an annulus for 
 * m=2 domains. It inherits from MultiplyConnectedDomain and provides the actual
 * boundary components, plus methods to update hole parameters during Newton iteration.
 *
 * Key features:
 * - Creates boundary components (unit circle + circular holes)
 * - Supports parameter updates during Newton iteration
 * - Provides initial guess computation from user domain analysis
 * - Handles both general case and annulus specialization
 */
class FornbergCanonicalDomain : public MultiplyConnectedDomain
{
public:
    /**
     * @brief Configuration type for the canonical domain
     */
    enum class CanonicalType
    {
        UNIT_DISK_WITH_HOLES,    // General m-connected: unit disk with (m-1) circular holes
        ANNULUS                  // Special 2-connected case: annulus between two circles
    };

    /**
     * @brief Initial guess strategy for conformal moduli
     */
    enum class InitialGuessStrategy
    {
        GEOMETRIC_CENTROIDS,     // Use boundary centroids and inscribed radii
        BOUNDARY_ANALYSIS,       // Analyze curvature and convex hull properties
        UNIFORM_DISTRIBUTION,    // Distribute holes uniformly in unit disk
        MANUAL                   // User-provided values
    };

private:
    // Canonical domain configuration
    CanonicalType m_canonical_type;
    std::vector<Complex> m_hole_centers;        // Centers c_ν of circular holes
    std::vector<double> m_hole_radii;           // Radii ρ_ν of circular holes
    
    // Boundary sampling
    int m_N;                                    // Points per boundary component

public:
    /**
     * @brief Construct canonical domain from hole parameters
     * @param hole_centers Centers c_ν of circular holes (empty for simply connected)
     * @param hole_radii Radii ρ_ν of circular holes (empty for simply connected)
     * @param N Number of boundary points per component (must be power of 2)
     */
    FornbergCanonicalDomain(const std::vector<Complex>& hole_centers = {},
                           const std::vector<double>& hole_radii = {},
                           int N = 128);

    /**
     * @brief Destructor
     */
    ~FornbergCanonicalDomain() = default;

    /**
     * @brief Update hole parameters during Newton iteration
     * @param new_centers New centers c_ν for circular holes
     * @param new_radii New radii ρ_ν for circular holes
     */
    void updateHoleParameters(const std::vector<Complex>& new_centers,
                             const std::vector<double>& new_radii);

    /**
     * @brief Get the canonical domain type
     * @return Canonical domain type (disk with holes vs. annulus)
     */
    CanonicalType getCanonicalType() const
    {
        return m_canonical_type;
    }

    /**
     * @brief Check if this is an annulus configuration
     * @return True if canonical domain is an annulus (m=2)
     */
    bool isAnnulus() const
    {
        return m_canonical_type == CanonicalType::ANNULUS;
    }

    /**
     * @brief Get hole centers c_ν
     * @return Vector of centers for circular holes
     */
    const std::vector<Complex>& getHoleCenters() const
    {
        return m_hole_centers;
    }

    /**
     * @brief Get hole radii ρ_ν
     * @return Vector of radii for circular holes
     */
    const std::vector<double>& getHoleRadii() const
    {
        return m_hole_radii;
    }

    /**
     * @brief Get conformal moduli as interleaved vector
     * @return Vector containing [c₁, ρ₁, c₂, ρ₂, ..., c_{m-1}, ρ_{m-1}]
     */
    Eigen::VectorXcd getConformalModuli() const;

    /**
     * @brief Set conformal moduli from interleaved vector
     * @param moduli Vector containing [c₁, ρ₁, c₂, ρ₂, ..., c_{m-1}, ρ_{m-1}]
     */
    void setConformalModuli(const Eigen::VectorXcd& moduli);

    /**
     * @brief Get number of boundary points per component
     * @return Number of points N
     */
    int getBoundaryPointCount() const
    {
        return m_N;
    }

    /**
     * @brief Validate that hole configuration is physically realizable
     * @return True if configuration is valid
     */
    bool isValidConfiguration() const;

    // Static factory methods for initial guess computation

    /**
     * @brief Compute initial guess from user domain geometric analysis
     * @param user_domain The user's multiply connected domain
     * @param strategy Strategy for computing initial guess
     * @param N Number of boundary points per component
     * @return New canonical domain with initial guess parameters
     */
    static std::shared_ptr<FornbergCanonicalDomain> 
        createFromUserDomain(std::shared_ptr<MultiplyConnectedDomain> user_domain,
                           InitialGuessStrategy strategy = InitialGuessStrategy::GEOMETRIC_CENTROIDS,
                           int N = 128);

    /**
     * @brief Compute geometric centroid-based initial guess
     * @param user_domain The user's multiply connected domain
     * @return Pair of (centers, radii) for initial guess
     */
    static std::pair<std::vector<Complex>, std::vector<double>>
        computeGeometricCentroidGuess(std::shared_ptr<MultiplyConnectedDomain> user_domain);

    /**
     * @brief Compute uniform distribution initial guess
     * @param connectivity Number of boundary components
     * @return Pair of (centers, radii) for initial guess
     */
    static std::pair<std::vector<Complex>, std::vector<double>>
        computeUniformDistributionGuess(int connectivity);

private:
    /**
     * @brief Create initial boundary components for constructor
     */
    static std::vector<std::shared_ptr<Boundary>> createInitialBoundaries(
        const std::vector<Complex>& hole_centers, 
        const std::vector<double>& hole_radii);

    /**
     * @brief Determine canonical type from connectivity
     */
    void determineCanonicalType();

    /**
     * @brief Create circular boundary component
     * @param center Center of the circle
     * @param radius Radius of the circle
     * @param is_outer True for outer boundary (counterclockwise), false for holes (clockwise)
     * @return Boundary component
     */
    static std::shared_ptr<Boundary> createCircularBoundary(const Complex& center, 
                                                           double radius, 
                                                           bool is_outer);
};