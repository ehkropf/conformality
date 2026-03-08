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
#include "../core/ConformalMap.h"
#include "../domains/Domain.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

MethodInfo TheodorsenMethod::getMethodInfo() const
{
    MethodInfo info;
    info.name = "Theodorsen";

    info.parameters = {
        {"Sample points", static_cast<int>(n_points)},
    };

    info.results = {
        {"Iterations", m_iteration_count},
        {"Residual", residual_norm},
        {"Converged", is_converged},
    };

    return info;
}

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

    // Store target domain for iteration (this is the starlike domain we map TO)
    m_current_domain = target_domain;

    // Validate domains - source should be unit circle, target should be starlike
    validateDomainCompatibility(source_domain, 1); // Check source connectivity
    validateSourceDomain(source_domain);           // Check source is unit circle
    validateDomainCompatibility(target_domain, 1); // Check target connectivity  
    validateTargetDomain(target_domain);           // Check target is starlike

    // Sample the boundary of the target domain (starlike)
    boundary_samples = sampleBoundary(target_domain, external);

    // Compute initial boundary moduli from target domain
    std::vector<double> rho = computeBoundaryModuli(target_domain, external);

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

void TheodorsenMethod::validateSourceDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("Source domain cannot be null");
    }

    auto circular_domain = std::dynamic_pointer_cast<CircularDomain>(domain);
    if (!circular_domain)
    {
        throw std::invalid_argument("Theodorsen's method requires the source domain to be a circular domain");
    }

    // Check if it's a unit circle centered at origin
    if (std::abs(circular_domain->getRadius() - 1.0) > 1e-12)
    {
        throw std::invalid_argument("Theodorsen's method requires the source domain to be a unit circle (radius = 1.0)");
    }

    if (std::abs(circular_domain->getCenter()) > 1e-12)
    {
        throw std::invalid_argument("Theodorsen's method requires the source domain to be centered at the origin");
    }
}

void TheodorsenMethod::validateTargetDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("Target domain cannot be null");
    }

    auto starlike_domain = std::dynamic_pointer_cast<StarlikeDomain>(domain);
    if (!starlike_domain)
    {
        throw std::invalid_argument("Theodorsen's method requires the target domain to be a starlike domain");
    }
}

Complex TheodorsenMethod::map(const Complex& z) const
{
    if (!is_converged)
    {
        throw std::runtime_error("Map has not been computed yet");
    }

    // Evaluate Laurent series for Theodorsen's method: g(z) = a_{-1}*z + a_0 + a_1/z + a_2/z^2 + ...
    // where g maps from unit disk U to starlike domain D
    // FFTW coefficients are in order: [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
    Complex result{0.0, 0.0};
    size_t N = laurent_coeffs.size();

    if (laurent_coeffs.empty())
    {
        return result;
    }

    // The DFT of boundary data ρ(φ)e^{iφ} gives us coefficients for Laurent series
    // For Theodorsen conformal map from unit disk to starlike domain:
    // g(z) = sum_{k=-∞}^{∞} c_k z^k where the map is analytic in the unit disk

    // Add constant term
    result += laurent_coeffs[0];

    // For mapping from unit disk to starlike domain, the dominant term should be z
    // Add z term (should be the dominant contribution, approximately = z for conformal map)
    if (N > 1) {
        // The coefficient for z is at index N-1 in FFTW ordering
        result += laurent_coeffs[N-1] * z;
    }

    // Add 1/z terms (Laurent series principal part)
    if (std::abs(z) > 1e-12) {
        Complex z_inv = Complex{1.0, 0.0} / z;
        Complex z_inv_power = z_inv;

        // Add first few 1/z^k terms only
        for (size_t k = 1; k < std::min(size_t(8), N/2); ++k) {
            result += laurent_coeffs[k] * z_inv_power;
            z_inv_power *= z_inv;
        }
    }

    return result;
}

Complex TheodorsenMethod::inverseMap([[maybe_unused]] const Complex& w) const
{
    if (!is_converged)
    {
        throw std::runtime_error("Map has not been computed yet");
    }
    
    // TODO: Implement inverse mapping from starlike domain back to unit disk
    // This requires either:
    // 1. Newton's method iteration to solve g(z) = w for z
    // 2. Barycentric Cauchy formula discretization (as mentioned in notebook)
    // 3. Using the boundary correspondence to interpolate back
    
    // For now, throw an exception with detailed information
    throw std::runtime_error("Inverse mapping from starlike domain to unit disk not yet implemented. "
                            "This requires iterative inversion of the Laurent series or "
                            "barycentric Cauchy formula implementation.");
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
    // Create log(ρ(φ)) for conjugation operator
    std::vector<Complex> log_rho(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        // Compute log(ρ(φ_k)) - real-valued since ρ > 0
        log_rho[i] = Complex(std::log(rho[i]), 0.0);
    }

    // Apply conjugation operator to log(ρ(φ_k))
    FFTWWrapper& fftw = FFTWWrapper::get_instance();
    std::vector<Complex> delta_fft = fftw.forward_fft(log_rho);

    // Apply conjugation operator: multiply by ±i for positive/negative frequencies
    size_t m = n_points / 2;
    for (size_t k = 1; k < m; ++k)
    {
        delta_fft[k] = delta_fft[k] * Complex(0.0, -1.0); // multiply by -i
    }
    for (size_t k = m; k < n_points; ++k)
    {
        delta_fft[k] = delta_fft[k] * Complex(0.0, 1.0);  // multiply by +i
    }

    // Get delta_{k+1} by inverse FFT (take real part)
    std::vector<Complex> delta_complex = fftw.backward_fft(delta_fft);
    std::vector<double> delta(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        delta[i] = delta_complex[i].real();
    }

    // Renormalize: φ_{k+1} = δ_{k+1} + θ - δ_{k+1}[0]
    double delta_0 = delta[0];
    std::vector<double> phi_new(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        phi_new[i] = delta[i] + phi_sequence[i] - delta_0;
    }

    // Compute new ρ values at updated φ positions
    auto target_domain = std::dynamic_pointer_cast<StarlikeDomain>(m_current_domain);
    std::vector<double> new_rho(n_points);
    for (size_t i = 0; i < n_points; ++i)
    {
        double radius = target_domain->getRadius(phi_new[i]);
        if (external)
        {
            new_rho[i] = 1.0 / radius;  // External mapping
        }
        else
        {
            new_rho[i] = radius;        // Internal mapping
        }
    }

    return new_rho;
}

void TheodorsenMethod::computeLaurentCoefficients(const std::vector<Complex>& boundary_data)
{
    // Use FFT to compute Fourier coefficients
    FFTWWrapper& fftw = FFTWWrapper::get_instance();
    std::vector<Complex> coeffs = fftw.forward_fft(boundary_data);

    // Store coefficients in Laurent series format with normalization
    laurent_coeffs.resize(n_points);
    double normalization = 1.0 / static_cast<double>(n_points);

    // Apply normalization to coefficients
    // FFTW returns coefficients in standard order: [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
    for (size_t i = 0; i < n_points; ++i)
    {
        laurent_coeffs[i] = coeffs[i] * normalization;
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
