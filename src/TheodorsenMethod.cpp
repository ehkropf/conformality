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

#include "TheodorsenMethod.h"
#include "ConformalMap.h"
#include "Domain.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>

TheodorsenMethod::TheodorsenMethod(size_t num_points)
    : ConformalMapMethod()
    , n_points(num_points)
    , is_converged(false)
    , residual_norm(0.0)
{
    if (!isPowerOfTwo(num_points))
    {
        throw std::invalid_argument("Number of points must be a power of 2 for FFT efficiency");
    }
    
    // Initialize phi sequence (equally spaced angles)
    phi_sequence.resize(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        phi_sequence[i] = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n_points);
    }
    
    m_max_iterations = 100;
}

void TheodorsenMethod::compute(ConformalMap& map_instance, double target_accuracy)
{
    auto source_domain = map_instance.getSourceDomain();
    auto target_domain = map_instance.getTargetDomain();
    MappingType mapping_type = map_instance.getMappingType();
    bool external = (mapping_type == MappingType::EXTERIOR_TO_INTERIOR);
    
    // Validate domains
    validateDomain(source_domain, 1); // Simply connected and starlike
    validateDomainCompatibility(target_domain, 1); // Target just needs connectivity check
    
    // Check if target domain is unit circle
    auto circular_target = std::dynamic_pointer_cast<CircularDomain>(target_domain);
    if (!circular_target || std::abs(circular_target->getRadius() - 1.0) > 1e-12 ||
        std::abs(circular_target->getCenter()) > 1e-12)
    {
        throw std::invalid_argument("Theodorsen's method requires unit circle as target domain");
    }
    
    // Sample the boundary
    boundary_samples = sampleBoundary(source_domain, external);
    
    // Compute initial boundary moduli
    std::vector<double> rho = computeBoundaryModuli(source_domain, external);
    
    // Iterative refinement
    m_iteration_count = 0;
    is_converged = false;
    residual_norm = std::numeric_limits<double>::infinity();
    
    while (m_iteration_count < m_max_iterations && !is_converged)
    {
        std::vector<double> old_rho = rho;
        rho = theodorsenIteration(rho, external);
        
        is_converged = checkConvergence(old_rho, rho, target_accuracy);
        ++m_iteration_count;
    }
    
    if (!is_converged)
    {
        throw std::runtime_error("Theodorsen's method failed to converge within maximum iterations");
    }
    
    // Compute final Laurent coefficients
    std::vector<Complex> boundary_correspondence(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        Complex exp_iphi(std::cos(phi_sequence[i]), std::sin(phi_sequence[i]));
        boundary_correspondence[i] = exp_iphi * Complex(rho[i], 0.0);
    }
    
    computeLaurentCoefficients(boundary_correspondence);
    m_achieved_accuracy = residual_norm;
}

void TheodorsenMethod::validateDomainGeometry(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("Domain cannot be null");
    }

    // Check if the domain is starlike
    auto starlike_domain = std::dynamic_pointer_cast<StarlikeDomain>(domain);
    if (!starlike_domain)
    {
        throw std::invalid_argument("Theodorsen's method requires a starlike domain");
    }
}

Complex TheodorsenMethod::map(const Complex& z) const
{
    if (!is_converged)
    {
        throw std::runtime_error("Map has not been computed yet");
    }
    
    // Evaluate Laurent series: f(z) = sum_{n=-∞}^{∞} a_n z^n
    // For numerical stability, we split into positive and negative powers
    Complex result{0.0, 0.0};
    
    // Handle n=0 coefficient separately
    if (!laurent_coeffs.empty())
    {
        result = result + laurent_coeffs[0];
    }
    
    // Positive powers: a_n z^n for n > 0
    Complex z_power = z;
    for (size_t n = 1; n < laurent_coeffs.size() / 2; ++n)
    {
        result = result + laurent_coeffs[n] * z_power;
        z_power = z_power * z;
    }
    
    // Negative powers: a_{-n} z^{-n} for n > 0
    if (std::abs(z) > 1e-12) // Avoid division by zero
    {
        Complex z_inv = Complex{1.0, 0.0} / z;
        Complex z_inv_power = z_inv;
        for (size_t n = 1; n < laurent_coeffs.size() / 2; ++n)
        {
            size_t neg_index = laurent_coeffs.size() / 2 + n;
            if (neg_index < laurent_coeffs.size())
            {
                result = result + laurent_coeffs[neg_index] * z_inv_power;
                z_inv_power = z_inv_power * z_inv;
            }
        }
    }
    
    return result;
}

Complex TheodorsenMethod::inverseMap([[maybe_unused]] const Complex& w) const
{
    // TODO: Implement inverse mapping using Newton's method or similar
    // For now, throw an exception as a placeholder
    throw std::runtime_error("Inverse mapping not yet implemented for Theodorsen's method");
}

void TheodorsenMethod::setNumPoints(size_t num_points)
{
    if (!isPowerOfTwo(num_points))
    {
        throw std::invalid_argument("Number of points must be a power of 2 for FFT efficiency");
    }
    
    n_points = num_points;
    
    // Reinitialize phi sequence
    phi_sequence.resize(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        phi_sequence[i] = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n_points);
    }
    
    // Clear previous computation results
    laurent_coeffs.clear();
    boundary_samples.clear();
    is_converged = false;
}

std::vector<Complex> TheodorsenMethod::sampleBoundary(std::shared_ptr<Domain> domain, [[maybe_unused]] bool external) const
{
    auto starlike_domain = std::dynamic_pointer_cast<StarlikeDomain>(domain);
    Complex center = starlike_domain->getCenter();
    
    std::vector<Complex> samples(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        double angle = phi_sequence[i];
        double radius = starlike_domain->getRadius(angle);
        
        Complex boundary_point = center + Complex(radius, 0.0) * Complex(std::cos(angle), std::sin(angle));
        samples[i] = boundary_point;
    }
    
    return samples;
}

std::vector<double> TheodorsenMethod::computeBoundaryModuli(std::shared_ptr<Domain> domain, bool external) const
{
    auto starlike_domain = std::dynamic_pointer_cast<StarlikeDomain>(domain);
    [[maybe_unused]] Complex center = starlike_domain->getCenter();
    
    std::vector<double> rho(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        double angle = phi_sequence[i];
        rho[i] = starlike_domain->getRadius(angle);
        
        if (external)
        {
            // For external mapping, invert the radius
            rho[i] = 1.0 / rho[i];
        }
    }
    
    return rho;
}

std::vector<double> TheodorsenMethod::theodorsenIteration(const std::vector<double>& rho, bool external)
{
    // Create complex boundary correspondence function
    std::vector<Complex> psi(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        Complex exp_iphi(std::cos(phi_sequence[i]), std::sin(phi_sequence[i]));
        psi[i] = Complex(rho[i], 0.0) * exp_iphi;
    }
    
    // Apply conjugation operator using FFTW wrapper
    FFTWWrapper& fftw = FFTWWrapper::get_instance();
    std::vector<Complex> K_psi = fftw.conjugation_operator(psi);
    
    // Update boundary correspondence
    std::vector<double> new_rho(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        Complex exp_iphi(std::cos(phi_sequence[i]), std::sin(phi_sequence[i]));
        
        if (external)
        {
            // External case: ψ_{n+1} = ψ_n - K[ψ_n]
            Complex updated_psi = psi[i] - K_psi[i];
            new_rho[i] = std::abs(updated_psi);
        }
        else
        {
            // Internal case: ψ_{n+1} = ψ_n + K[ψ_n]
            Complex updated_psi = psi[i] + K_psi[i];
            new_rho[i] = std::abs(updated_psi);
        }
    }
    
    return new_rho;
}

void TheodorsenMethod::computeLaurentCoefficients(const std::vector<Complex>& boundary_data)
{
    // Use FFT to compute Fourier coefficients
    FFTWWrapper& fftw = FFTWWrapper::get_instance();
    std::vector<Complex> coeffs = fftw.forward_fft(boundary_data);
    
    // Store coefficients in Laurent series format
    laurent_coeffs.resize(n_points);
    
    // Copy coefficients with proper indexing
    // FFTW returns coefficients in standard order: [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
    for (size_t i = 0; i < n_points; ++i)
    {
        laurent_coeffs[i] = coeffs[i];
    }
}

bool TheodorsenMethod::checkConvergence(const std::vector<double>& old_rho, const std::vector<double>& new_rho, double tolerance)
{
    // Compute L2 norm of the difference
    double sum_sq_diff = 0.0;
    double sum_sq_old = 0.0;
    
    for (size_t i = 0; i < n_points; ++i)
    {
        double diff = new_rho[i] - old_rho[i];
        sum_sq_diff += diff * diff;
        sum_sq_old += old_rho[i] * old_rho[i];
    }
    
    residual_norm = std::sqrt(sum_sq_diff) / std::sqrt(sum_sq_old);
    return residual_norm < tolerance;
}

bool TheodorsenMethod::isPowerOfTwo(size_t num_points) const
{
    return num_points > 0 && (num_points & (num_points - 1)) == 0;
}