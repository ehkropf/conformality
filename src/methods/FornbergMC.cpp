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

#include "FornbergMC.h"
#include "PMatrixBuilder.h"
#include "../numerics/CGSolver.h"
#include "../numerics/Polyval.h"
#include "../domains/FornbergCanonicalDomain.h"
#include "../core/ConformalMap.h"
#include "../domains/Domain.h"
// #include "../core/StatusManager.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

FornbergMC::FornbergMC(const FornbergMCConfiguration& config)
    : ConformalMapMethod()
    , m_config{config}
{
    // Validate configuration
    m_config.validate();
    
    // Validate that N is a power of 2 for FFT compatibility
    if (!isPowerOfTwo(m_config.N))
    {
        throw std::invalid_argument("FornbergMC: Number of boundary points N must be a power of 2");
    }
}

void FornbergMC::compute(ConformalMap& map_instance, double target_accuracy)
{
    // Use the more restrictive tolerance between target_accuracy and configured newton_tolerance
    double effective_tolerance = std::min(target_accuracy, m_config.newton_tolerance);
    
    // TODO: Log computation start with domain connectivity and target accuracy
    
    // Validate domains
    validateDomains(map_instance);
    
    // Extract target domain
    auto target_domain = map_instance.getTargetDomain();
    mp_user_domain = std::dynamic_pointer_cast<MultiplyConnectedDomain>(target_domain);
    
    if (!mp_user_domain)
    {
        throw std::invalid_argument("FornbergMC: Target domain must be multiply connected");
    }
    
    m_connectivity = mp_user_domain->getConnectivity();
    
    // Create canonical domain from user domain analysis
    mp_canonical_domain = FornbergCanonicalDomain::createFromUserDomain(
        mp_user_domain,
        FornbergCanonicalDomain::InitialGuessStrategy::GEOMETRIC_CENTROIDS,
        m_config.N
    );
    
    // Detect annulus case and initialize
    m_is_annulus = detectAnnulusCase();
    if (m_is_annulus && m_config.verbose)
    {
        // TODO: Log "Detected annulus case (m=2) - using optimized algorithm"
    }
    
    // Initialize Newton iteration framework
    initializeNewtonIteration();
    
    // Newton iteration loop
    m_is_converged = false;
    m_residual_history.clear();
    
    for (int iter = 0; iter < m_config.max_newton_iterations && !m_is_converged; ++iter)
    {
        if (m_config.verbose)
        {
            printIterationDiagnostics(iter);
        }
        
        formSystem();
        solveSystem();
        newtonUpdate();
        m_is_converged = checkConvergence(effective_tolerance);
        
        // Store residual for history
        m_residual_history.push_back(m_current_residual);
        
        // Boundary redistribution if needed
        if (iter % m_config.redistribution_frequency == 0 && m_config.enable_redistribution)
        {
            if (redistributeBoundaryParameters())
            {
                // TODO: Log "Redistributed boundary parameters at iteration" + iter
            }
        }
    }
    
    if (!m_is_converged)
    {
        // TODO: Log warning "Newton iteration failed to converge after" + iterations + "iterations"
        if (!m_config.enable_fallback_methods)
        {
            throw std::runtime_error("FornbergMC: Newton iteration failed to converge");
        }
    }
    else
    {
        // TODO: Log "Newton iteration converged in" + iterations + "iterations, residual =" + m_current_residual
    }
    
    computeFourierCoefficients();
    
    m_achieved_accuracy = m_current_residual;
}

Complex FornbergMC::map(const Complex& z) const
{
    // Evaluate series representation matching MATLAB bmap.map_eval and bdd_eval
    // w = polyval(a(N/2:-1:1,1),z) + sum_nu polyval([a(N/2+1:N,nu);0], rho/(z-c))
    
    if (m_a.size() == 0)
    {
        throw std::runtime_error("FornbergMC: Map not computed yet");
    }
    
    if (m_connectivity == 0)
    {
        throw std::runtime_error("FornbergMC: Invalid connectivity state");
    }
    
    const int N = m_a.rows();
    const int m = m_a.cols();
    
    if (N % 2 != 0)
    {
        throw std::runtime_error("FornbergMC: N must be even for coefficient matrix splitting");
    }
    
    const int half_N = N / 2;
    
    // First part: Taylor series using polyval(a(N/2:-1:1,1), z)
    // Extract coefficients a(N/2:-1:1,1) in descending order for polyval
    std::vector<Complex> taylor_coeffs;
    taylor_coeffs.reserve(half_N);
    
    for (int j = half_N - 1; j >= 0; --j) // N/2-1, N/2-2, ..., 0 (0-based)
    {
        taylor_coeffs.push_back(m_a(j, 0));
    }
    
    Complex result = taylor_coeffs.empty() ? Complex(0.0, 0.0) : polyval(taylor_coeffs, z);
    
    // Second part: Laurent series for each hole nu = 2 to m  
    // polyval([a(N/2+1:N,nu);0], rho(nu-1)/(z-c(nu-1)))
    if (m_connectivity > 1 && mp_canonical_domain)
    {
        const auto& hole_centers = mp_canonical_domain->getHoleCenters();
        const auto& hole_radii = mp_canonical_domain->getHoleRadii();
        
        for (int nu = 1; nu < m; ++nu) // nu = 2,3,...,m in MATLAB (1,2,...,m-1 in 0-based)
        {
            if (nu - 1 >= static_cast<int>(hole_centers.size())) break;
            
            const Complex c_nu = hole_centers[nu - 1];
            const double rho_nu = hole_radii[nu - 1];
            
            const Complex z_minus_c = z - c_nu;
            
            // Check for singularity at hole center
            if (std::abs(z_minus_c) < 1e-14)
            {
                throw std::runtime_error("FornbergMC: Evaluation point too close to hole center");
            }
            
            const Complex w_arg = rho_nu / z_minus_c; // Argument to polyval: rho/(z-c)
            
            // MATLAB: [a(N/2+1:N,nu);0] creates coefficients with appended zero
            // This is [a(half_N, nu), a(half_N+1, nu), ..., a(N-1, nu), 0]
            std::vector<Complex> laurent_coeffs;
            laurent_coeffs.reserve(half_N + 1);
            
            for (int j = half_N; j < N; ++j) // half_N, half_N+1, ..., N-1
            {
                laurent_coeffs.push_back(m_a(j, nu));
            }
            laurent_coeffs.push_back(Complex(0.0, 0.0)); // Append zero as in MATLAB
            
            result += laurent_coeffs.empty() ? Complex(0.0, 0.0) : polyval(laurent_coeffs, w_arg);
        }
    }
    
    return result;
}

Complex FornbergMC::inverseMap(const Complex& w) const
{
    // Stub implementation - requires iterative solution
    if (m_a.size() == 0)
    {
        throw std::runtime_error("FornbergMC: Map not computed yet");
    }
    
    // TODO: Implement inverse map evaluation (typically via Newton's method)
    // TODO: Log warning "Inverse map evaluation not implemented - returning identity"
    return w; // Placeholder
}

void FornbergMC::setConfiguration(const FornbergMCConfiguration& config)
{
    config.validate();
    m_config = config;
    m_is_converged = false;
    // TODO: Log configuration update with key parameters (N, tolerances, etc.)
}

void FornbergMC::validateSourceDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("FornbergMC: Source domain cannot be null");
    }
    
    // Source domain should be canonical (unit disk with circular holes)
    if (domain->getConnectivity() < 2)
    {
        throw std::invalid_argument("FornbergMC: Source domain must be multiply connected (connectivity >= 2)");
    }
    
    // TODO: Add more specific validation for canonical domain structure?
}

void FornbergMC::validateTargetDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("FornbergMC: Target domain cannot be null");
    }
    
    if (domain->getConnectivity() < 2)
    {
        throw std::invalid_argument("FornbergMC: Target domain must be multiply connected (connectivity >= 2)");
    }
    
    if (domain->isUnbounded())
    {
        throw std::invalid_argument("FornbergMC: Target domain must be bounded");
    }
}

void FornbergMC::initializeNewtonIteration()
{
    // TODO: Log "Initializing Newton iteration with N=" + m_config.N + ", connectivity=" + m_connectivity
    
    mp_matrix_builder = std::make_unique<PMatrixBuilder>(m_config, m_connectivity, m_is_annulus);
    mp_cg_solver = std::make_unique<CGSolver>(m_config);
    
    sampleBoundaries();
    
    initializeConformalModuli();
    
    // Initialize matrices and vectors
    int system_size = m_connectivity * m_config.N; // Approximate size
    m_S.resize(m_connectivity, m_config.N);
    m_conformal_moduli.resize(2 * (m_connectivity - 1)); // c_ν, ρ_ν pairs
    m_D.resize(system_size, system_size);
    m_g.resize(system_size);
    m_U.resize(system_size);
    m_a.resize(m_config.N, m_connectivity);
    
    // TODO: Log "Initialized system matrices: D(" + system_size + "x" + system_size + "), coefficients a(" + m_config.N + "x" + m_connectivity + ")"
}

void FornbergMC::formSystem()
{
    // Stub implementation for system formation
    // TODO: Debug log "Forming linear system D*U = g"
    
    // TODO: Implement system formation using PMatrixBuilder
    // This involves:
    // 1. Constructing P_ν matrices for analyticity conditions
    // 2. Setting up discretized boundary conditions
    // 3. Handling annulus case optimization if applicable
    // 4. Forming D matrix and g vector
}

void FornbergMC::solveSystem()
{
    // Stub implementation for system solution
    // TODO: Debug log "Solving linear system with CG method"
    
    // TODO: Implement CG solution using CGSolver
    // This involves:
    // 1. Setting up function handle for D†D*x = D†*g
    // 2. Calling custom CG solver with convergence monitoring
    // 3. Handling best-iterate tracking for non-convergent cases
}

void FornbergMC::newtonUpdate()
{
    // TODO: Implement Newton update of boundary correspondences and conformal moduli
    // This involves:
    // 1. Extracting updates for S_ν from solution vector U
    // 2. Updating conformal moduli c_ν, ρ_ν
    // 3. Applying damping if enabled
    // 4. Computing residual norm for convergence check
    
    // For now, just update the canonical domain with current moduli
    if (mp_canonical_domain && m_conformal_moduli.size() > 0)
    {
        try
        {
            mp_canonical_domain->setConformalModuli(m_conformal_moduli);
        }
        catch (const std::exception&)
        {
            // If update fails, keep current parameters
        }
    }
    
    // Placeholder: set some residual value
    m_current_residual = 1e-6; // This should be computed from actual updates
}

bool FornbergMC::checkConvergence(double tolerance)
{
    // Check Newton update norm convergence
    bool newton_converged = m_current_residual < tolerance;
    if (m_config.verbose && newton_converged)
    {
        // TODO: Log "Newton iteration converged: residual = " + m_current_residual
    }
    return newton_converged;
}

bool FornbergMC::detectAnnulusCase() const
{
    if (!m_config.auto_detect_annulus)
    {
        return m_config.force_annulus_mode;
    }
    return (m_connectivity == 2);
}

void FornbergMC::initializeConformalModuli()
{
    if (!mp_canonical_domain)
    {
        throw std::runtime_error("FornbergMC: Canonical domain must be created before initializing conformal moduli");
    }
    
    m_conformal_moduli = mp_canonical_domain->getConformalModuli();
}

void FornbergMC::sampleBoundaries()
{
    // TODO: Debug log "Sampling " + m_connectivity + " boundary components with N=" + m_config.N + " points each"
    
    m_boundary_samples.resize(m_connectivity);
    m_parameter_values.resize(m_connectivity);
    
    // TODO: Sample each boundary component at N parameter values
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        m_boundary_samples[nu].resize(m_config.N);
        m_parameter_values[nu].resize(m_config.N);
        
        // Placeholder: uniform parameter sampling
        for (int j = 0; j < m_config.N; ++j)
        {
            double theta = 2.0 * M_PI * j / m_config.N;
            m_parameter_values[nu][j] = theta;
            m_boundary_samples[nu][j] = Complex(std::cos(theta), std::sin(theta)); // Unit circle placeholder
        }
    }
}

bool FornbergMC::redistributeBoundaryParameters()
{
    // Stub implementation for boundary parameter redistribution
    // TODO: Debug log "Checking boundary parameter distribution quality"
    
    // TODO: Implement redistribution logic
    // This involves:
    // 1. Checking for parameter clustering or gaps
    // 2. Redistributing parameters if quality metrics exceeded
    // 3. Resampling boundaries with new parameter values
    
    return false; // No redistribution performed
}

void FornbergMC::computeFourierCoefficients()
{
    // TODO: Debug log "Computing Fourier coefficients via FFT"
    
    // TODO: Implement Fourier coefficient computation from boundary data
    // This involves FFTs of boundary correspondence functions
}

bool FornbergMC::isPowerOfTwo(int N) const
{
    return N > 0 && (N & (N - 1)) == 0;
}

void FornbergMC::printIterationDiagnostics(size_t iteration) const
{
    // TODO: Log iteration diagnostics: "Iteration " + iteration + ": residual = " + m_current_residual
    
    if (m_config.eigenvalue_analysis)
    {
        // TODO: Print eigenvalue diagnostics
    }
}