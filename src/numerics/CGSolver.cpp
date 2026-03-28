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

#include "CGSolver.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

CGSolver::CGSolver(const FornbergMCConfiguration& config)
    : m_config{config}
    , m_best_residual{std::numeric_limits<double>::max()}
    , m_best_iteration{-1}
    , m_restart_count{0}
{
}

Eigen::VectorXd CGSolver::solve(const MatrixVectorProduct& A_function, 
                               const Eigen::VectorXd& b, 
                               const Eigen::VectorXd& x0) const
{
    // Validate input
    validateSystem(A_function, b, x0);
    
    // Initialize solver state
    initializeSolverState(b.size());
    
    // Determine initial guess
    Eigen::VectorXd initial_guess = x0.size() > 0 ? x0 : Eigen::VectorXd::Zero(b.size());
    
    // Run CG iteration
    ConvergenceInfo info;
    Eigen::VectorXd solution = cgIteration(A_function, b, initial_guess, info, m_config.max_cgm_iterations);
    
    // Store convergence information
    info.restart_count = m_restart_count;
    m_last_convergence_info = info;

    // Log final convergence status
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    if (info.converged)
    {
        oss << "Converged in " << info.iterations << " iterations, "
            << "relative residual = " << info.relative_residual;
        mp_statusManager->reportInfo("CGSolver", oss.str());
    }
    else
    {
        oss << "Did not converge after " << info.iterations << " iterations, "
            << "relative residual = " << info.relative_residual;
        if (info.used_best_iterate)
        {
            oss << " (using best iterate from iteration " << info.best_iterate_index << ")";
        }
        mp_statusManager->reportWarning("CGSolver", oss.str());
    }

    return solution;
}

Eigen::VectorXd CGSolver::solveComplexSystem(
    const std::function<Eigen::VectorXcd(const Eigen::VectorXd&)>& D_function,
    const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& D_adjoint_function,
    const Eigen::VectorXcd& g,
    const Eigen::VectorXd& x0) const
{
    // Create real system function: 2*real(D†D)*x
    auto real_A_function = createRealSystemFunction(D_function, D_adjoint_function);

    // Convert complex RHS to real: 2*real(g)
    Eigen::VectorXd real_b = 2.0 * g.real();

    // Solve the real system
    return solve(real_A_function, real_b, x0);
}

bool CGSolver::runSelfTest(int size) const
{
    mp_statusManager->reportDebug("CGSolver", "Running self-test with system size " + std::to_string(size));

    try
    {
        // Create a simple test system: A = I + 0.1*ones, x_true = ones, b = A*x_true
        Eigen::VectorXd x_true = Eigen::VectorXd::Ones(size);

        auto test_A_function = [size](const Eigen::VectorXd& x) -> Eigen::VectorXd
        {
            return x + 0.1 * x.sum() * Eigen::VectorXd::Ones(size);
        };

        Eigen::VectorXd b = test_A_function(x_true);
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(size);

        // Solve the test system
        Eigen::VectorXd x_computed = solve(test_A_function, b, x0);

        // Check error
        double error = (x_computed - x_true).norm();
        bool passed = error < 1e-10;

        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6);
        oss << "Self-test " << (passed ? "passed" : "failed") << ", error = " << error;
        if (passed)
        {
            mp_statusManager->reportInfo("CGSolver", oss.str());
        }
        else
        {
            mp_statusManager->reportWarning("CGSolver", oss.str());
        }

        return passed;
    }
    catch (const std::exception& e)
    {
        mp_statusManager->reportError("CGSolver", "Self-test failed with exception", e.what());
        return false;
    }
}

Eigen::VectorXd CGSolver::cgIteration(const MatrixVectorProduct& A_function,
                                     const Eigen::VectorXd& b,
                                     const Eigen::VectorXd& x0,
                                     ConvergenceInfo& info,
                                     int max_iterations) const
{
    const int n = b.size();
    const double b_norm = b.norm();
    
    if (b_norm == 0.0)
    {
        // Trivial case: zero RHS
        info.converged = true;
        info.iterations = 0;
        info.final_residual = 0.0;
        info.relative_residual = 0.0;
        return Eigen::VectorXd::Zero(n);
    }

    // Early exit on non-finite RHS (avoids running all CG iterations with NaN/inf propagation)
    if (!std::isfinite(b_norm))
    {
        std::string msg = "Non-finite RHS norm detected (" + std::to_string(b_norm) + "). Aborting CG.";
        mp_statusManager->reportWarning("CGSolver", msg);
        info.converged = false;
        info.iterations = 0;
        info.final_residual = b_norm;
        info.relative_residual = std::numeric_limits<double>::infinity();
        return x0;
    }

    // Initialize CG variables
    Eigen::VectorXd x = x0;
    Eigen::VectorXd r = b - A_function(x);
    Eigen::VectorXd p = r;

    double rsold = r.dot(r);
    double initial_residual = std::sqrt(rsold);

    // Early exit on non-finite initial residual
    if (!std::isfinite(initial_residual))
    {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6);
        oss << "Non-finite initial residual (" << initial_residual << "). Aborting CG.";
        mp_statusManager->reportWarning("CGSolver", oss.str());
        info.converged = false;
        info.iterations = 0;
        info.final_residual = initial_residual;
        info.relative_residual = std::numeric_limits<double>::infinity();
        return x0;
    }

    info.residual_history.clear();
    info.residual_history.reserve(max_iterations);
    info.residual_history.push_back(initial_residual);

    // Initialize best iterate tracking (MATLAB cgm.m line 19: x_ = r)
    if (m_config.enable_best_iterate)
    {
        updateBestIterate(r, initial_residual, 0);
    }

    // Log initial state
    {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6);
        oss << "Starting CG: n=" << n << ", initial residual=" << initial_residual
            << ", b_norm=" << b_norm;
        mp_statusManager->reportDebug("CGSolver", oss.str());
    }

    // CG iteration loop
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // Compute A*p
        Eigen::VectorXd Ap = A_function(p);
        double pAp = p.dot(Ap);

        // Safety-net: catch degenerate systems (negative-definite, NaN/Inf).
        // MATLAB cgm.m has no breakdown check — for well-posed Fornberg systems,
        // D'D is positive definite and pAp > 0 throughout. This guard catches
        // truly degenerate or corrupted systems.
        if (pAp <= 0.0 || !std::isfinite(pAp))
        {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(6);
            oss << "CG breakdown at iteration " << iter << ": p'Ap = " << pAp
                << " (expected positive for SPD system)";
            mp_statusManager->reportWarning("CGSolver", oss.str());
            info.converged = false;
            info.iterations = iter;
            info.final_residual = std::sqrt(rsold);
            info.relative_residual = info.final_residual / b_norm;
            if (m_config.enable_best_iterate && m_best_iteration >= 0)
            {
                info.used_best_iterate = true;
                info.best_iterate_index = m_best_iteration;
                return m_best_iterate;
            }
            info.used_best_iterate = false;
            return x;
        }

        // Update solution and residual
        double alpha = rsold / pAp;
        x = x + alpha * p;
        r = r - alpha * Ap;
        
        double rsnew = r.dot(r);
        double current_residual = std::sqrt(rsnew);
        double relative_residual = current_residual / b_norm;

        // Early exit on non-finite residual (NaN/inf cascade)
        if (!std::isfinite(current_residual))
        {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(6);
            oss << "Non-finite residual at iteration " << iter << ": " << current_residual;
            mp_statusManager->reportWarning("CGSolver", oss.str());
            info.converged = false;
            info.iterations = iter + 1;
            info.final_residual = current_residual;
            info.relative_residual = std::numeric_limits<double>::infinity();
            if (m_config.enable_best_iterate && m_best_iteration >= 0)
            {
                info.used_best_iterate = true;
                info.best_iterate_index = m_best_iteration;
                return m_best_iterate;
            }
            info.used_best_iterate = false;
            return x;
        }

        info.residual_history.push_back(current_residual);
        // Update best iterate
        if (m_config.enable_best_iterate)
        {
            updateBestIterate(x, current_residual, iter + 1);
        }
        
        // Check convergence
        if (relative_residual < m_config.cgm_tolerance)
        {
            info.converged = true;
            info.iterations = iter + 1;
            info.final_residual = current_residual;
            info.relative_residual = relative_residual;
            info.used_best_iterate = false;
            return x;
        }
        
        // Check for restart condition
        if (shouldRestart(info.residual_history, iter + 1))
        {
            {
                std::ostringstream oss;
                oss << std::scientific << std::setprecision(6);
                oss << "CG restart triggered at iteration " << (iter + 1)
                    << ", residual=" << current_residual;
                mp_statusManager->reportDebug("CGSolver", oss.str());
            }
            info.iterations = iter + 1;
            return performRestart(A_function, b, x, max_iterations - (iter + 1), info);
        }
        
        // Update search direction
        double beta = rsnew / rsold;
        p = r + beta * p;
        rsold = rsnew;
    }
    
    // Did not converge within iteration limit
    info.converged = false;
    info.iterations = max_iterations;
    info.final_residual = std::sqrt(rsold);
    info.relative_residual = info.final_residual / b_norm;
    
    // Return best iterate if enabled and available
    if (m_config.enable_best_iterate && m_best_iteration >= 0)
    {
        info.used_best_iterate = true;
        info.best_iterate_index = m_best_iteration;
        return m_best_iterate;
    }
    else
    {
        info.used_best_iterate = false;
        return x;
    }
}

bool CGSolver::shouldRestart(const std::vector<double>& residual_history, int current_iter) const
{
    if (m_config.max_cgm_restarts == 0)
    {
        return false; // Restarts disabled
    }

    if (current_iter < 10 || residual_history.size() < 10)
    {
        return false; // Too early to judge stagnation
    }
    
    // Check if residual has stagnated
    const int check_length = 5;
    if (residual_history.size() < check_length + 1)
    {
        return false;
    }
    
    double recent_residual = residual_history.back();
    double old_residual = residual_history[residual_history.size() - check_length - 1];
    
    double improvement_rate = (old_residual - recent_residual) / old_residual;
    
    return improvement_rate < m_config.cgm_restart_threshold;
}

Eigen::VectorXd CGSolver::performRestart(const MatrixVectorProduct& A_function,
                                        const Eigen::VectorXd& b,
                                        const Eigen::VectorXd& current_x,
                                        int remaining_iters,
                                        ConvergenceInfo& info) const
{
    if (remaining_iters <= 0)
    {
        return current_x;
    }

    if (m_restart_count >= m_config.max_cgm_restarts)
    {
        std::string msg = "Maximum restart count (" + std::to_string(m_config.max_cgm_restarts) + ") reached";
        mp_statusManager->reportWarning("CGSolver", msg);
        return current_x;
    }
    m_restart_count++;

    // Save current convergence info
    auto saved_info = info;
    
    // Restart CG with current iterate as initial guess
    ConvergenceInfo restart_info;
    Eigen::VectorXd restarted_solution = cgIteration(A_function, b, current_x, restart_info, remaining_iters);
    
    // Merge convergence information
    info.iterations = saved_info.iterations + restart_info.iterations;
    info.converged = restart_info.converged;
    info.final_residual = restart_info.final_residual;
    info.relative_residual = restart_info.relative_residual;
    
    // Merge residual histories
    info.residual_history.insert(info.residual_history.end(),
                                restart_info.residual_history.begin(),
                                restart_info.residual_history.end());
    
    return restarted_solution;
}

void CGSolver::updateBestIterate(const Eigen::VectorXd& x, double residual, int iteration) const
{
    if (residual < m_best_residual)
    {
        double improvement = (m_best_iteration >= 0) ? (m_best_residual - residual) / m_best_residual : 1.0;
        m_best_residual = residual;
        m_best_iterate = x;
        m_best_iteration = iteration;

        // Log significant improvements in best iterate (when enabled and significant)
        if (iteration == 0 || improvement > 0.1)
        {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(6);
            oss << "New best iterate at iteration " << iteration << ", residual=" << residual;
            mp_statusManager->reportDebug("CGSolver", oss.str());
        }
    }
}

void CGSolver::initializeSolverState(int system_size) const
{
    m_best_residual = std::numeric_limits<double>::max();
    m_best_iteration = -1;
    m_best_iterate.resize(system_size);
    m_restart_count = 0;
}

void CGSolver::validateSystem(const MatrixVectorProduct& A_function,
                             const Eigen::VectorXd& b,
                             const Eigen::VectorXd& x0) const
{
    if (b.size() == 0)
    {
        throw std::invalid_argument("CGSolver: Empty RHS vector");
    }
    
    if (x0.size() > 0 && x0.size() != b.size())
    {
        throw std::invalid_argument("CGSolver: Initial guess size mismatch");
    }
    
    // Test that A_function works with the given vector size
    try
    {
        Eigen::VectorXd test_x = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd test_result = A_function(test_x);
        
        if (test_result.size() != b.size())
        {
            throw std::invalid_argument("CGSolver: Matrix-vector product size mismatch");
        }
    }
    catch (const std::exception& e)
    {
        throw std::invalid_argument("CGSolver: Matrix-vector product function failed: " + std::string(e.what()));
    }
}

CGSolver::MatrixVectorProduct CGSolver::createRealSystemFunction(
    const std::function<Eigen::VectorXcd(const Eigen::VectorXd&)>& D_function,
    const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& D_adjoint_function) const
{
    return [D_function, D_adjoint_function](const Eigen::VectorXd& x) -> Eigen::VectorXd
    {
        // Compute D*x (real input, complex output)
        Eigen::VectorXcd Dx = D_function(x);
        
        // Compute D†(D*x) (complex input, complex output)  
        Eigen::VectorXcd DtDx = D_adjoint_function(Dx);
        
        // Return 2*real(D†D*x) as real vector
        return 2.0 * DtDx.real();
    };
}
