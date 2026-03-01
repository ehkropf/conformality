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

#include "../core/StatusManager.h"
#include "../methods/FornbergMCConfiguration.h"
#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <vector>

/**
 * @brief Custom conjugate gradient solver for the Fornberg method linear systems
 *
 * This solver is specifically designed for the Fornberg-like method's linear systems
 * of the form D†D*x = D†g, where D is the system matrix from analyticity conditions.
 * The key innovation is using a function handle approach to compute D†D*x as D†(D*x)
 * without ever forming the dense matrix D†D explicitly.
 *
 * Key features:
 * - Function handle approach avoids O(N²) memory usage
 * - Best iterate tracking for non-convergent cases
 * - Complex arithmetic support via real system transformation
 * - Custom restart strategies and convergence monitoring
 * - Direct port from MATLAB fornmc/cgm.m implementation
 */
class CGSolver
{
public:
    /**
     * @brief Function type for matrix-vector products
     * Takes a vector and returns the result of applying a linear operator
     */
    using MatrixVectorProduct = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

    /**
     * @brief Convergence information structure
     */
    struct ConvergenceInfo
    {
        bool converged = false;
        int iterations = 0;
        double final_residual = 0.0;
        double relative_residual = 0.0;
        std::vector<double> residual_history;
        bool used_best_iterate = false;
        int best_iterate_index = -1;
        int restart_count = 0;
    };

private:
    const FornbergMCConfiguration& m_config;
    std::shared_ptr<IStatusManager> mp_statusManager;

    // Solver state
    mutable ConvergenceInfo m_last_convergence_info;
    mutable Eigen::VectorXd m_best_iterate;
    mutable double m_best_residual;
    mutable int m_best_iteration;
    mutable int m_restart_count;

public:
    /**
     * @brief Construct a new CG Solver
     * @param config Configuration parameters
     */
    explicit CGSolver(const FornbergMCConfiguration& config);

    /**
     * @brief Destructor
     */
    ~CGSolver() = default;

    /**
     * @brief Solve the linear system A*x = b using conjugate gradient method
     * @param A_function Function handle for matrix-vector product A*x
     * @param b Right-hand side vector
     * @param x0 Initial guess (optional, zero if not provided)
     * @return Solution vector x
     */
    Eigen::VectorXd solve(const MatrixVectorProduct& A_function, 
                         const Eigen::VectorXd& b, 
                         const Eigen::VectorXd& x0 = Eigen::VectorXd()) const;

    /**
     * @brief Solve the complex system 2*real(D†D)*x = 2*real(D†g)
     * This is the main interface for Fornberg method linear systems
     * @param D_function Function handle for D*x (complex matrix-vector product)
     * @param D_adjoint_function Function handle for D†*x (adjoint matrix-vector product)
     * @param g Complex right-hand side vector D†*g_original
     * @param x0 Initial guess for real solution vector
     * @return Real solution vector x
     */
    Eigen::VectorXd solveComplexSystem(const std::function<Eigen::VectorXcd(const Eigen::VectorXd&)>& D_function,
                                      const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& D_adjoint_function,
                                      const Eigen::VectorXcd& g,
                                      const Eigen::VectorXd& x0 = Eigen::VectorXd()) const;

    /**
     * @brief Get convergence information from the last solve
     * @return Convergence information structure
     */
    const ConvergenceInfo& getLastConvergenceInfo() const
    {
        return m_last_convergence_info;
    }

    /**
     * @brief Check if the last solve converged
     * @return True if converged
     */
    bool hasConverged() const
    {
        return m_last_convergence_info.converged;
    }

    /**
     * @brief Get the final residual from the last solve
     * @return Final residual norm
     */
    double getFinalResidual() const
    {
        return m_last_convergence_info.final_residual;
    }

    /**
     * @brief Get the number of iterations from the last solve
     * @return Iteration count
     */
    int getIterationCount() const
    {
        return m_last_convergence_info.iterations;
    }

    /**
     * @brief Test the solver with a simple known system
     * @param size System size for test
     * @return True if test passed
     */
    bool runSelfTest(int size = 100) const;

    /**
     * @brief Set the status manager for logging
     * @param manager Shared pointer to status manager (nullptr to disable logging)
     */
    void setStatusManager(std::shared_ptr<IStatusManager> manager)
    {
        mp_statusManager = manager;
    }

    /**
     * @brief Get the current status manager
     * @return Shared pointer to the status manager (may be null)
     */
    std::shared_ptr<IStatusManager> getStatusManager() const
    {
        return mp_statusManager;
    }

private:
    /**
     * @brief Core CG iteration implementation
     * @param A_function Matrix-vector product function
     * @param b Right-hand side vector
     * @param x0 Initial guess
     * @param info Convergence information (output)
     * @param max_iterations Maximum number of CG iterations for this invocation
     * @return Solution vector
     */
    Eigen::VectorXd cgIteration(const MatrixVectorProduct& A_function,
                               const Eigen::VectorXd& b,
                               const Eigen::VectorXd& x0,
                               ConvergenceInfo& info,
                               int max_iterations) const;

    /**
     * @brief Check if restart is needed based on stagnation
     * @param residual_history History of residual norms
     * @param current_iter Current iteration
     * @return True if restart should be performed
     */
    bool shouldRestart(const std::vector<double>& residual_history, int current_iter) const;

    /**
     * @brief Perform CG restart with current iterate as new starting point
     *
     * Returns current_x without restarting if the maximum restart count has been
     * reached or no iteration budget remains.
     *
     * @param A_function Matrix-vector product function
     * @param b Right-hand side vector
     * @param current_x Current iterate
     * @param remaining_iters Remaining iteration budget
     * @param info Convergence information (input/output)
     * @return Restarted solution vector, or current_x if restart was skipped
     */
    Eigen::VectorXd performRestart(const MatrixVectorProduct& A_function,
                                  const Eigen::VectorXd& b,
                                  const Eigen::VectorXd& current_x,
                                  int remaining_iters,
                                  ConvergenceInfo& info) const;

    /**
     * @brief Update best iterate tracking
     * @param x Current iterate
     * @param residual Current residual norm
     * @param iteration Current iteration number
     */
    void updateBestIterate(const Eigen::VectorXd& x, double residual, int iteration) const;

    /**
     * @brief Initialize solver state for a new solve
     * @param system_size Size of the system
     */
    void initializeSolverState(int system_size) const;

    /**
     * @brief Validate system compatibility
     * @param A_function Matrix-vector product function
     * @param b Right-hand side vector
     * @param x0 Initial guess
     */
    void validateSystem(const MatrixVectorProduct& A_function,
                       const Eigen::VectorXd& b,
                       const Eigen::VectorXd& x0) const;

    /**
     * @brief Convert complex system to real form for CG solution
     * D†D is Hermitian positive definite, so 2*real(D†D) is real symmetric positive definite
     * @param D_function Function for D*x
     * @param D_adjoint_function Function for D†*x  
     * @return Function handle for 2*real(D†D)*x
     */
    MatrixVectorProduct createRealSystemFunction(
        const std::function<Eigen::VectorXcd(const Eigen::VectorXd&)>& D_function,
        const std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>& D_adjoint_function) const;
};
