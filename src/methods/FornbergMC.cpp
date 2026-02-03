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
#include "../numerics/FFTWWrapper.h"
#include "../domains/FornbergCanonicalDomain.h"
#include "../core/ConformalMap.h"
#include "../domains/Domain.h"
// #include "../core/StatusManager.h"
#include <iostream>
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

    // Correct system dimensions (D matrix is rectangular, not square)
    const int N = m_config.N;
    const int M = N / 2;
    const int m = m_connectivity;

    int num_rows, num_cols;
    if (m_is_annulus)
    {
        num_rows = m * M;
        num_cols = m * N + 1;  // Only rho for annulus (c fixed at 0)
    }
    else
    {
        num_rows = m * M + 2;              // +2 normalization rows
        num_cols = m * N + 3 * (m - 1);    // +rho, +Re(c), +Im(c) per inner boundary
    }

    m_S.resize(N, m);  // (N, m) layout to match MATLAB column-per-boundary convention
    m_conformal_moduli.resize(2 * (m - 1));
    m_D.resize(num_rows, num_cols);
    m_g.resize(num_rows);
    m_U.resize(num_cols);
    m_abs_eta.resize(N, m);  // (N, m) layout to match MATLAB
    m_a.resize(N, m);

    // Initialize S to identity (uniform theta spacing)
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            m_S(j, nu) = 2.0 * M_PI * j / N;
        }
    }

    // TODO: Log "Initialized system matrices: D(" + num_rows + "x" + num_cols + "), coefficients a(" + N + "x" + m + ")"
}

void FornbergMC::formSystem()
{
    // Validate preconditions
    if (!mp_user_domain)
    {
        throw std::runtime_error(
            "FornbergMC::formSystem: User domain not initialized. Call compute() first.");
    }
    if (!mp_matrix_builder)
    {
        throw std::runtime_error(
            "FornbergMC::formSystem: Matrix builder not initialized. Call initializeNewtonIteration() first.");
    }

    const int N = m_config.N;
    const int M = N / 2;
    const int m = m_connectivity;

    // Validate boundary count matches connectivity
    const auto& boundaries = mp_user_domain->getBoundaries();
    if (static_cast<int>(boundaries.size()) != m)
    {
        throw std::runtime_error(
            "FornbergMC::formSystem: Boundary count mismatch. Expected " + std::to_string(m) +
            " boundaries, got " + std::to_string(boundaries.size()) + ".");
    }

    // Zero out system
    m_D.setZero();
    m_g.setZero();

    // Fourier grid: q_k = exp(-i * 2*pi*k / N)
    Eigen::VectorXcd q(N);
    for (int k = 0; k < N; ++k)
    {
        double theta = -2.0 * M_PI * k / N;
        q(k) = Complex(std::cos(theta), std::sin(theta));
    }

    // Sample target boundary at current S values
    Eigen::MatrixXcd xi(N, m);   // boundary positions
    Eigen::MatrixXcd eta(N, m);  // normalized tangents

    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            double S_j = m_S(j, nu);
            xi(j, nu) = boundaries[nu]->evaluate(S_j);
            Complex tangent = boundaries[nu]->evaluateDerivative(S_j);
            double abs_t = std::abs(tangent);
            m_abs_eta(j, nu) = abs_t;
            if (abs_t > 1e-14)
            {
                eta(j, nu) = tangent / abs_t;
            }
            else
            {
                // TODO: Replace with spdlog warning when logging is integrated
                std::cerr << "Warning: Degenerate tangent detected at boundary " << nu
                          << ", point " << j << " (S=" << S_j << "). Using fallback unit tangent."
                          << std::endl;
                eta(j, nu) = Complex(1.0, 0.0);
            }
        }
    }

    // Build ConformalModuli struct for PMatrixBuilder
    ConformalModuli moduli;
    moduli.c.resize(m - 1);
    moduli.rho.resize(m - 1);
    for (int i = 0; i < m - 1; ++i)
    {
        moduli.c(i) = m_conformal_moduli(2 * i);
        moduli.rho(i) = std::real(m_conformal_moduli(2 * i + 1));
    }

    // Build normalization rows p1/pnu
    // p1: constructed here matching MATLAB namap.m constructor (lines 49-52)
    // pnu: computed per inner boundary (MATLAB namap.m form_system, lines 92-93)
    // Only needed for general (non-annulus) case
    Eigen::RowVectorXcd p1;
    std::vector<Eigen::RowVectorXcd> pnu_vectors;

    if (!m_is_annulus)
    {
        // Using z_0 = 0 (normalization point at origin)
        constexpr Complex z_0(0.0, 0.0);

        // p1 for outer boundary: [1, z_0, z_0^2, ..., z_0^(M-1), 0, 0, ..., 0]
        // When z_0 = 0: p1 = [1, 0, 0, ..., 0]
        // Note: Non-zero z_0 support would require generalizing this initialization
        p1 = Eigen::RowVectorXcd::Zero(N);
        p1(0) = 1.0;

        // pnu vectors for inner boundaries (one per inner boundary)
        // pnu[M:N-1] = (rho/(z_0-c))^(M:-1:1)
        pnu_vectors.resize(m - 1);
        for (int nu_idx = 0; nu_idx < m - 1; ++nu_idx)
        {
            pnu_vectors[nu_idx] = Eigen::RowVectorXcd::Zero(N);
            Complex c_val = moduli.c(nu_idx);
            double rho_val = moduli.rho(nu_idx);

            // Validate that inner boundary center is not at normalization point
            Complex denom = z_0 - c_val;
            if (std::abs(denom) < 1e-14)
            {
                throw std::invalid_argument(
                    "FornbergMC::formSystem: Inner boundary center c(" + std::to_string(nu_idx) +
                    ") is too close to normalization point z_0. Domain may be degenerate.");
            }
            Complex ratio = rho_val / denom;  // = rho / (-c) when z_0 = 0

            for (int k = 0; k < M; ++k)
            {
                // MATLAB: pnu(M+1:N) = ratio.^(M:-1:1)
                // C++ 0-based: pnu(M+k) = ratio^(M-k) for k=0..M-1
                pnu_vectors[nu_idx](M + k) = std::pow(ratio, M - k);
            }
        }
    }

    FFTWWrapper& fftw = FFTWWrapper::get_instance();

    // Number of rows populated by P_nu blocks (before final normalization row is added):
    // Annulus: m*M rows (no p1/pnu extension, no additional normalization rows)
    // General: m*M+1 rows (P_nu extended with p1/pnu; applyNormalization adds the final row to reach m*M+2)
    const int num_constraint_rows = m_is_annulus ? (m * M) : (m * M + 1);

    for (int nu = 0; nu < m; ++nu)
    {
        // Get base P_nu from builder (m*M × N)
        Eigen::MatrixXcd P_base = mp_matrix_builder->buildPMatrix(nu, moduli);

        Eigen::MatrixXcd P_nu;
        if (m_is_annulus)
        {
            // Annulus case: use P_base directly (no p1/pnu extension)
            P_nu = P_base;
        }
        else
        {
            // General case: extend with p1 (outer) or pnu (inner) to get (m*M+1 × N)
            P_nu.resize(num_constraint_rows, N);
            P_nu.topRows(m * M) = P_base;
            if (nu == 0)
            {
                P_nu.row(m * M) = p1;
            }
            else
            {
                P_nu.row(m * M) = pnu_vectors[nu - 1];
            }
        }

        int col_offset = nu * N;

        // Build D block: P_nu * FFT(diag(eta_nu))
        // Each column k is P_nu * FFT(e_k * eta(k,nu))
        for (int k = 0; k < N; ++k)
        {
            std::vector<Complex> col_k(N, Complex(0.0, 0.0));
            col_k[k] = eta(k, nu);
            std::vector<Complex> fft_col = fftw.forward_fft(col_k);
            Eigen::VectorXcd fft_vec = Eigen::Map<Eigen::VectorXcd>(fft_col.data(), N);
            m_D.block(0, col_offset + k, num_constraint_rows, 1) = P_nu * fft_vec;
        }

        // RHS contribution: g -= P_nu * FFT(xi_nu)
        std::vector<Complex> xi_col(N);
        for (int j = 0; j < N; ++j) xi_col[j] = xi(j, nu);
        std::vector<Complex> fft_xi = fftw.forward_fft(xi_col);
        Eigen::VectorXcd fft_xi_vec = Eigen::Map<Eigen::VectorXcd>(fft_xi.data(), N);
        m_g.head(num_constraint_rows) -= P_nu * fft_xi_vec;

        // For inner boundaries, add moduli derivative columns
        if (nu >= 1)
        {
            double rho_nu = moduli.rho(nu - 1);

            // Validate rho is positive (non-degenerate inner boundary)
            if (rho_nu <= 0.0)
            {
                throw std::runtime_error(
                    "FornbergMC::formSystem: Invalid conformal modulus rho(" + std::to_string(nu - 1) +
                    ") = " + std::to_string(rho_nu) + ". Inner boundary radius must be positive.");
            }

            // Compute S derivative: Sdiff(j) = (S(j+1) - S(j)) * N/(2*pi)
            Eigen::VectorXd S_diff(N);
            for (int j = 0; j < N; ++j)
            {
                int j_next = (j + 1) % N;
                S_diff(j) = (m_S(j_next, nu) - m_S(j, nu)) * N / (2.0 * M_PI);
            }

            // zeta = i * |eta| * eta * Sdiff / rho
            std::vector<Complex> zeta(N);
            for (int j = 0; j < N; ++j)
            {
                zeta[j] = Complex(0, 1) * m_abs_eta(j, nu) * eta(j, nu) * S_diff(j) / rho_nu;
            }
            std::vector<Complex> fft_zeta = fftw.forward_fft(zeta);
            Eigen::VectorXcd fft_zeta_vec = Eigen::Map<Eigen::VectorXcd>(fft_zeta.data(), N);

            // Column for rho derivative
            int rho_col = m * N + (nu - 1);
            m_D.block(0, rho_col, num_constraint_rows, 1) = P_nu * fft_zeta_vec;

            if (!m_is_annulus)
            {
                // Columns for c derivative (Re and Im parts)
                std::vector<Complex> q_zeta(N);
                for (int j = 0; j < N; ++j) q_zeta[j] = q(j) * zeta[j];
                std::vector<Complex> fft_q_zeta = fftw.forward_fft(q_zeta);
                Eigen::VectorXcd fft_q_zeta_vec = Eigen::Map<Eigen::VectorXcd>(fft_q_zeta.data(), N);

                int re_c_col = m * N + (m - 1) + 2 * (nu - 1);
                int im_c_col = re_c_col + 1;
                m_D.block(0, re_c_col, num_constraint_rows, 1) = P_nu * fft_q_zeta_vec;
                m_D.block(0, im_c_col, num_constraint_rows, 1) = Complex(0, 1) * P_nu * fft_q_zeta_vec;
            }
        }
    }

    // Apply normalization (adds constraint rows for general case)
    // Using default norm_cond = [1, 0, 0] from MATLAB reference:
    // - z_0 = 0 (normalization point at origin, affects p1/pnu via constructor)
    // - norm_value = 0 (added to g(end-1))
    if (!m_is_annulus)
    {
        constexpr double norm_value = 0.0;  // norm_cond(3) = 0
        mp_matrix_builder->applyNormalizationConditions(m_D, m_g, norm_value);
    }
}

void FornbergMC::solveSystem()
{
    // Precondition checks
    if (!mp_cg_solver)
    {
        throw std::runtime_error(
            "FornbergMC::solveSystem: CGSolver not initialized");
    }

    if (m_D.rows() == 0 || m_g.size() == 0)
    {
        throw std::runtime_error(
            "FornbergMC::solveSystem: System not formed");
    }

    // Create D_function: applies D to a real vector
    auto D_function = [this](const Eigen::VectorXd& x) -> Eigen::VectorXcd {
        return m_D * x;
    };

    // Create D_adjoint_function: applies D† to a complex vector
    auto D_adjoint_function = [this](const Eigen::VectorXcd& y) -> Eigen::VectorXcd {
        return m_D.adjoint() * y;
    };

    // Compute transformed RHS: D† * g
    Eigen::VectorXcd g_transformed = m_D.adjoint() * m_g;

    // Solve the system
    Eigen::VectorXd solution = mp_cg_solver->solveComplexSystem(
        D_function,
        D_adjoint_function,
        g_transformed
    );

    // Store solution in m_U (as complex with zero imaginary part)
    m_U = solution.cast<std::complex<double>>();

    // Log convergence information
    const auto& info = mp_cg_solver->getLastConvergenceInfo();

    if (!info.converged)
    {
        std::cerr << "Warning: CG solver did not converge after "
                  << info.iterations << " iterations. "
                  << "Residual: " << info.final_residual;
        if (info.used_best_iterate)
        {
            std::cerr << " (using best iterate from iteration "
                      << info.best_iterate_index << ")";
        }
        std::cerr << '\n';
    }
    else if (m_config.verbose)
    {
        std::cout << "CG converged in " << info.iterations
                  << " iterations, residual: " << info.final_residual
                  << '\n';
    }
}

void FornbergMC::newtonUpdate()
{
    // Compute residual as infinity norm of update vector (matches MATLAB implementation)
    // This must be done BEFORE applying updates to state variables
    m_current_residual = m_U.lpNorm<Eigen::Infinity>();

    const int m = m_connectivity;
    const int N = m_config.N;

    // Step 1: Scale boundary updates by 1/abs_eta (matches MATLAB: U(1:m*N)./abs_eta(:))
    constexpr double kEpsilon = 1e-14;
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            double abs_eta_val = m_abs_eta(j, nu);
            if (abs_eta_val > kEpsilon)
            {
                m_U(nu * N + j) /= abs_eta_val;
            }
            else
            {
                // TODO: Replace with spdlog warning when logging is integrated
                std::cerr << "Warning: Degenerate abs_eta detected at boundary " << nu
                          << ", point " << j << " (abs_eta=" << abs_eta_val
                          << "). Skipping scaling for this point." << std::endl;
            }
        }
    }

    // Step 2: Determine damping factor
    double damping = m_config.enable_newton_damping ? m_config.newton_damping_factor : 1.0;

    // Step 3: Update boundary correspondences S_ν(θ)
    // MATLAB: obj.S = obj.S + reshape(obj.U(1:m*N), N, m)
    for (int nu = 0; nu < m; ++nu)
    {
        for (int j = 0; j < N; ++j)
        {
            m_S(j, nu) += damping * std::real(m_U(nu * N + j));
        }
    }

    // Step 4: Update conformal moduli
    // Solution vector layout after S updates:
    //   m*N to m*N+m-2: radii updates (δρ₂, δρ₃, ..., δρ_m)
    //   m*N+m-1 onwards: center updates as interleaved Re/Im pairs

    // Update radii (all cases): solution[m*N+k] → moduli[2*k+1]
    // MATLAB: obj.rho = obj.rho + obj.U(m*N+1:m*N+m-1)
    for (int k = 0; k < m - 1; ++k)
    {
        int sol_idx = m * N + k;
        int mod_idx = 2 * k + 1;  // Radii are at odd indices in interleaved storage
        double radius_update = damping * std::real(m_U(sol_idx));
        m_conformal_moduli(mod_idx) += radius_update;

        // Validate radius remains positive
        if (std::real(m_conformal_moduli(mod_idx)) <= 0.0)
        {
            throw std::runtime_error(
                "FornbergMC::newtonUpdate: Radius for boundary " + std::to_string(k + 2) +
                " became non-positive (" + std::to_string(std::real(m_conformal_moduli(mod_idx))) +
                "). Newton iteration has diverged - consider enabling damping or checking domain configuration.");
        }
    }

    // Update centers (m >= 3 only): solution[m*N+m-1+2k, m*N+m+2k] → moduli[2*k]
    // MATLAB: obj.c = obj.c + obj.U(m*N+m:2:end) + 1i*obj.U(m*N+m+1:2:end)
    if (!m_is_annulus)
    {
        for (int k = 0; k < m - 1; ++k)
        {
            int re_idx = m * N + m - 1 + 2 * k;  // Real part of δc_{k+2}
            int im_idx = m * N + m + 2 * k;      // Imaginary part of δc_{k+2}
            int mod_idx = 2 * k;                  // Centers are at even indices

            Complex center_update(
                damping * std::real(m_U(re_idx)),
                damping * std::real(m_U(im_idx))
            );
            m_conformal_moduli(mod_idx) += center_update;
        }
    }
    // Annulus case (m=2): c₂ stays at 0, no center updates needed

    // Step 5: Sync canonical domain with updated moduli
    if (mp_canonical_domain && m_conformal_moduli.size() > 0)
    {
        try
        {
            mp_canonical_domain->setConformalModuli(m_conformal_moduli);
        }
        catch (const std::runtime_error& e)
        {
            // Only catch runtime errors - validation errors (std::invalid_argument)
            // should propagate as they indicate programming bugs
            // TODO: Replace with spdlog::warn() when logging infrastructure is integrated
            std::cerr << "Warning: Failed to update canonical domain - "
                      << e.what() << ". Keeping current parameters." << std::endl;
        }
    }
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

    if (!mp_user_domain)
    {
        throw std::runtime_error("FornbergMC: User domain must be set before sampling boundaries");
    }

    const auto& boundaries = mp_user_domain->getBoundaries();

    m_boundary_samples.resize(m_connectivity);
    m_parameter_values.resize(m_connectivity);

    // Sample each boundary component at N uniformly spaced parameter values
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        m_boundary_samples[nu].resize(m_config.N);
        m_parameter_values[nu].resize(m_config.N);

        for (int j = 0; j < m_config.N; ++j)
        {
            double theta = 2.0 * M_PI * j / m_config.N;
            m_parameter_values[nu][j] = theta;
            m_boundary_samples[nu][j] = boundaries[nu]->evaluate(theta, 0);
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

    // Get FFT instance and normalization factor
    FFTWWrapper& fftw = FFTWWrapper::get_instance();
    const double normalization = 1.0 / static_cast<double>(m_config.N);

    // Get user domain boundaries
    const auto& boundaries = mp_user_domain->getBoundaries();

    // For each boundary component, compute Fourier coefficients from boundary evaluations
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        // Step 1: Evaluate user boundary at converged S_ν parameters
        // S_ν(θ_j) gives the parameter on user boundary ν corresponding to canonical parameter θ_j
        std::vector<Complex> xi(m_config.N);
        for (int j = 0; j < m_config.N; ++j)
        {
            // Get boundary correspondence parameter
            double S_nu_j = m_S(j, nu);

            // Evaluate user boundary at this parameter to get the mapped point
            xi[j] = boundaries[nu]->evaluate(S_nu_j, 0);
        }

        // Step 2: Compute FFT of boundary points
        std::vector<Complex> coeffs = fftw.forward_fft(xi);

        // Step 3: Normalize and store in m_a matrix (column nu)
        // FFTW returns unnormalized coefficients - divide by N for proper DFT
        for (int j = 0; j < m_config.N; ++j)
        {
            m_a(j, nu) = coeffs[j] * normalization;
        }
    }
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
