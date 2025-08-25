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

#include "FornbergCanonicalDomain.h"
#include "BoundaryComponent.h"
#include "Boundary.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

FornbergCanonicalDomain::FornbergCanonicalDomain(const std::vector<Complex>& hole_centers,
                                                 const std::vector<double>& hole_radii,
                                                 int N)
    : MultiplyConnectedDomain(createInitialBoundaries(hole_centers, hole_radii), false)
    , m_hole_centers{hole_centers}
    , m_hole_radii{hole_radii}
    , m_N{N}
{
    // Validate input
    if (hole_centers.size() != hole_radii.size())
    {
        throw std::invalid_argument("FornbergCanonicalDomain: Number of centers must match number of radii");
    }
    
    if (N <= 0 || (N & (N - 1)) != 0)
    {
        throw std::invalid_argument("FornbergCanonicalDomain: N must be a positive power of 2");
    }
    
    // Validate hole parameters
    for (double r : hole_radii)
    {
        if (r <= 0.0)
        {
            throw std::invalid_argument("FornbergCanonicalDomain: All hole radii must be positive");
        }
    }
    
    // Determine canonical type
    determineCanonicalType();
    
    // Validate geometric configuration
    if (!isValidConfiguration())
    {
        throw std::invalid_argument("FornbergCanonicalDomain: Invalid hole configuration - holes overlap or extend outside unit disk");
    }
}

void FornbergCanonicalDomain::updateHoleParameters(const std::vector<Complex>& new_centers,
                                                   const std::vector<double>& new_radii)
{
    if (new_centers.size() != m_hole_centers.size() || new_radii.size() != m_hole_radii.size())
    {
        throw std::invalid_argument("FornbergCanonicalDomain: Cannot change number of holes during update");
    }
    
    // Validate new parameters
    for (double r : new_radii)
    {
        if (r <= 0.0)
        {
            throw std::invalid_argument("FornbergCanonicalDomain: All hole radii must be positive");
        }
    }
    
    // Update parameters
    m_hole_centers = new_centers;
    m_hole_radii = new_radii;
    
    // Recreate boundary components with new parameters
    boundaries = createInitialBoundaries(m_hole_centers, m_hole_radii);
}

Eigen::VectorXcd FornbergCanonicalDomain::getConformalModuli() const
{
    int num_holes = static_cast<int>(m_hole_centers.size());
    Eigen::VectorXcd moduli(2 * num_holes);
    
    for (int i = 0; i < num_holes; ++i)
    {
        moduli[2*i] = m_hole_centers[i];        // Center c_ν
        moduli[2*i + 1] = Complex(m_hole_radii[i], 0.0);  // Radius ρ_ν (real)
    }
    
    return moduli;
}

void FornbergCanonicalDomain::setConformalModuli(const Eigen::VectorXcd& moduli)
{
    int num_holes = static_cast<int>(m_hole_centers.size());
    if (moduli.size() != 2 * num_holes)
    {
        throw std::invalid_argument("FornbergCanonicalDomain: Moduli size must be 2 * number of holes");
    }
    
    std::vector<Complex> new_centers(num_holes);
    std::vector<double> new_radii(num_holes);
    
    for (int i = 0; i < num_holes; ++i)
    {
        new_centers[i] = moduli[2*i];
        new_radii[i] = std::abs(moduli[2*i + 1]);  // Take magnitude to ensure positive radius
        
        if (new_radii[i] <= 0.0)
        {
            throw std::invalid_argument("FornbergCanonicalDomain: All radii must be positive");
        }
    }
    
    updateHoleParameters(new_centers, new_radii);
}

bool FornbergCanonicalDomain::isValidConfiguration() const
{
    // Check that radii are positive
    for (double r : m_hole_radii)
    {
        if (r <= 0.0)
        {
            return false;
        }
    }
    
    // Check that holes don't overlap with unit circle or each other
    for (size_t i = 0; i < m_hole_centers.size(); ++i)
    {
        double r_i = m_hole_radii[i];
        Complex c_i = m_hole_centers[i];
        
        // Check distance from origin (unit circle boundary)
        if (std::abs(c_i) + r_i >= 1.0)
        {
            return false; // Hole extends outside unit disk
        }
        
        // Check overlap with other holes
        for (size_t j = i + 1; j < m_hole_centers.size(); ++j)
        {
            double r_j = m_hole_radii[j];
            Complex c_j = m_hole_centers[j];
            
            if (std::abs(c_i - c_j) < r_i + r_j)
            {
                return false; // Holes overlap
            }
        }
    }
    
    return true;
}

std::shared_ptr<FornbergCanonicalDomain> 
FornbergCanonicalDomain::createFromUserDomain(std::shared_ptr<MultiplyConnectedDomain> user_domain,
                                              InitialGuessStrategy strategy,
                                              int N)
{
    if (!user_domain)
    {
        throw std::invalid_argument("FornbergCanonicalDomain: User domain cannot be null");
    }
    
    int connectivity = user_domain->getConnectivity();
    if (connectivity < 2)
    {
        throw std::invalid_argument("FornbergCanonicalDomain: User domain must be multiply connected");
    }
    
    std::vector<Complex> centers;
    std::vector<double> radii;
    
    switch (strategy)
    {
        case InitialGuessStrategy::GEOMETRIC_CENTROIDS:
            std::tie(centers, radii) = computeGeometricCentroidGuess(user_domain);
            break;
            
        case InitialGuessStrategy::UNIFORM_DISTRIBUTION:
            std::tie(centers, radii) = computeUniformDistributionGuess(connectivity);
            break;
            
        case InitialGuessStrategy::BOUNDARY_ANALYSIS:
            // FIXME: For now, fall back to geometric centroids
            std::tie(centers, radii) = computeGeometricCentroidGuess(user_domain);
            break;
            
        case InitialGuessStrategy::MANUAL:
            throw std::invalid_argument("FornbergCanonicalDomain: Manual strategy requires explicit parameters");
    }
    
    return std::make_shared<FornbergCanonicalDomain>(centers, radii, N);
}

std::pair<std::vector<Complex>, std::vector<double>>
FornbergCanonicalDomain::computeGeometricCentroidGuess(std::shared_ptr<MultiplyConnectedDomain> user_domain)
{
    int connectivity = user_domain->getConnectivity();
    int num_holes = connectivity - 1;
    
    std::vector<Complex> centers(num_holes);
    std::vector<double> radii(num_holes);
    
    // Simple heuristic: distribute holes based on boundary component analysis
    const auto& boundaries = user_domain->getBoundaries();
    
    for (int i = 0; i < num_holes; ++i)
    {
        // For interior boundaries (holes), place them at scaled centroids
        int boundary_idx = i + 1; // Skip outer boundary at index 0
        
        if (boundary_idx < static_cast<int>(boundaries.size()))
        {
            // FIXME: this is a placeholder
            // Compute a simple centroid estimate
            // In practice, this would analyze the actual boundary geometry
            double angle = 2.0 * M_PI * i / num_holes;
            double scale_factor = 0.3; // Conservative placement inside unit disk
            
            centers[i] = scale_factor * Complex(std::cos(angle), std::sin(angle));
            radii[i] = 0.1; // Conservative radius
        }
        else
        {
            // Fallback to uniform distribution
            double angle = 2.0 * M_PI * i / num_holes;
            centers[i] = 0.3 * Complex(std::cos(angle), std::sin(angle));
            radii[i] = 0.1;
        }
    }
    
    return std::make_pair(centers, radii);
}

std::pair<std::vector<Complex>, std::vector<double>>
FornbergCanonicalDomain::computeUniformDistributionGuess(int connectivity)
{
    int num_holes = connectivity - 1;
    
    std::vector<Complex> centers(num_holes);
    std::vector<double> radii(num_holes);
    
    // Distribute holes uniformly in a circle inside unit disk
    double distribution_radius = 0.4; // Stay well inside unit circle
    double uniform_hole_radius = 0.08; // Small uniform holes
    
    for (int i = 0; i < num_holes; ++i)
    {
        double angle = 2.0 * M_PI * i / num_holes;
        centers[i] = distribution_radius * Complex(std::cos(angle), std::sin(angle));
        radii[i] = uniform_hole_radius;
    }
    
    return std::make_pair(centers, radii);
}


void FornbergCanonicalDomain::determineCanonicalType()
{
    int num_holes = static_cast<int>(m_hole_centers.size());
    
    if (num_holes == 0)
    {
        throw std::invalid_argument("FornbergCanonicalDomain: Cannot create simply connected canonical domain");
    }
    else if (num_holes == 1)
    {
        m_canonical_type = CanonicalType::ANNULUS;
    }
    else
    {
        m_canonical_type = CanonicalType::UNIT_DISK_WITH_HOLES;
    }
}

std::vector<std::shared_ptr<Boundary>> FornbergCanonicalDomain::createInitialBoundaries(
    const std::vector<Complex>& hole_centers,
    const std::vector<double>& hole_radii)
{
    std::vector<std::shared_ptr<Boundary>> boundaries;
    
    // Create outer boundary (unit circle, counterclockwise)
    boundaries.push_back(createCircularBoundary(Complex(0.0, 0.0), 1.0, true));
    
    // Create hole boundaries (clockwise orientation)
    for (size_t i = 0; i < hole_centers.size(); ++i)
    {
        boundaries.push_back(createCircularBoundary(hole_centers[i], hole_radii[i], false));
    }
    
    return boundaries;
}

std::shared_ptr<Boundary> FornbergCanonicalDomain::createCircularBoundary(const Complex& center,
                                                                         double radius,
                                                                         bool is_outer)
{
    // Create an analytic boundary component for a circle
    auto boundary_component = std::make_shared<AnalyticBoundaryComponent>(
        [center, radius, is_outer](double theta) -> Complex {
            // Parameterization: z(θ) = center + radius * e^(iθ)
            // For holes (inner boundaries), reverse orientation with -θ
            double actual_theta = is_outer ? theta : -theta;
            return center + radius * Complex(std::cos(actual_theta), std::sin(actual_theta));
        },
        [radius, is_outer](double theta) -> Complex {
            // Derivative: z'(θ) = radius * i * e^(iθ) * (±1)
            double actual_theta = is_outer ? theta : -theta;
            double sign = is_outer ? 1.0 : -1.0;
            return sign * radius * Complex(-std::sin(actual_theta), std::cos(actual_theta));
        }
    );
    
    return std::make_shared<Boundary>(std::vector<std::shared_ptr<BoundaryComponent>>{boundary_component});
}