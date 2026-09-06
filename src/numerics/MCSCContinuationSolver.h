/*
 * Copyright © 2026, Everett Kropf (ehkropf@gmail.com)
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

#include <Eigen/Dense>
#include <functional>
#include <stdexcept>

/**
 * @brief Numerical continuation solver for F(X) = 0 (dissertation Sec 3.3; port of numcontin.m's
 *        continuation(), Algower & Georg's Contup / Program 3).
 *
 * Not MCSC-specific: operates on any smooth residual callable F: R^n -> R^n. Solves F(X) = 0 by
 * tracing the zero curve of the homotopy H(X, lambda) = F(X) + (lambda - 1)*F(X0) from lambda = 0
 * (trivially satisfied at X = X0) to lambda = 1 (F(X) = 0), via an Euler predictor / chord-Newton
 * corrector method with a Broyden-updated finite-difference Jacobian and adaptive step-size
 * control.
 *
 * Chosen (per project design decision, #160) over adapting FornbergMC's damped-Newton approach,
 * because arbitrary MCSC polygon/circle-domain pairs can be geometrically far apart and
 * continuation is markedly more robust to a poor initial guess -- the dissertation (Sec 3.3.2,
 * Table 3.1) documents MATLAB's fsolve/lsqnonlin/fminunc/fminsearch all failing or badly
 * underperforming continuation on exactly this parameter problem.
 *
 * Deliberately drops numcontin.m's monitor/plot/figure-handle machinery (no GUI story yet) and its
 * fprintf progress output (no stdout side effects from a library class).
 */
class MCSCContinuationSolver
{
public:
    using ResidualFunction = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

    /**
     * @brief Tunable continuation parameters, mirroring numcontin.m's property defaults.
     */
    struct Options
    {
        /// Target infinity-norm of the homotopy residual. numcontin.m defaults to 1e-15, but that
        /// is tighter than the accumulated floating-point noise from the forward-difference
        /// Jacobian (finiteDifferenceStep ~ 1e-6) and repeated Broyden updates can reliably reach
        /// once lambda snaps to 1 in the final Newton phase -- 1e-13 is the achievable analog.
        double tolerance = 1e-13;
        double hmin = 1e-14;           ///< Minimum step size before giving up.
        double hmax = 1.0;             ///< Maximum step size.
        double resmin = 1e-8;          ///< Floor added to the corrector-contraction denominator.
        double contractionMax = 0.5;   ///< Corrector contraction ratio above which a step is halved.
        double angleMax = 0.78539816339744830961;  ///< pi/4: max tangent-turning angle per step.
        double initialStep = 1.0 / 16.0;
        int maxRefinementsBeforeJacobianReset = 3;
        double finiteDifferenceStep = 1e-6;
        int maxSteps = 100000;  ///< Safety cap on predictor-corrector steps (not in numcontin.m).
    };

    /**
     * @brief Result of a continuation solve.
     */
    struct Result
    {
        Eigen::VectorXd solution;
        double residualNorm = 0.0;
        double lastLambda = 0.0;
        double lastStepSize = 0.0;
    };

    /**
     * @brief Exception thrown when continuation cannot reach lambda = 1 within tolerance.
     *
     * Raised when the step size underflows Options::hmin before the homotopy residual reaches
     * Options::tolerance, or when Options::maxSteps predictor-corrector steps are exhausted --
     * both cases numcontin.m leaves as a silent early exit from its while loop, which this
     * project's error-handling convention (CLAUDE.md: solver loops must not silently swallow a
     * non-convergent result) instead signals as a hard failure.
     */
    class ConvergenceError : public std::runtime_error
    {
    public:
        explicit ConvergenceError(const std::string& message)
            : std::runtime_error(message)
        {
        }
    };

    /**
     * @throws std::invalid_argument if any Options field is out of its valid range (see
     *         validateOptions()).
     */
    explicit MCSCContinuationSolver(ResidualFunction F);
    MCSCContinuationSolver(ResidualFunction F, Options options);

    /**
     * @brief Run continuation from x0 to a zero of F.
     * @param x0 Initial guess (F(x0) need not be near zero).
     * @return Solution and diagnostic information.
     * @throws ConvergenceError if the solve does not reach Options::tolerance before the step
     *         size underflows Options::hmin or Options::maxSteps steps are exhausted, if the
     *         homotopy residual or Jacobian becomes non-finite at any step (e.g. F evaluated a
     *         degenerate input), or if the continuation Jacobian becomes rank-deficient (a
     *         near-zero diagonal entry in its QR factorization) so no reliable tangent direction
     *         exists.
     */
    Result solve(const Eigen::VectorXd& x0) const;

private:
    ResidualFunction m_F;
    Options m_options;

    /// @throws std::invalid_argument if any field of options is out of its valid range (a
    /// non-positive tolerance/step-size/count, or hmax <= hmin).
    static void validateOptions(const Options& options);

    /// Homotopy H([lambda; X]) = F(X) + (lambda - 1)*F0, matching numcontin.m's curvefun.
    Eigen::VectorXd curveFunction(const Eigen::VectorXd& extendedState, const Eigen::VectorXd& F0) const;

    /// Forward-difference Jacobian of curveFunction at extendedState (n x (n+1)), port of numjac.
    Eigen::MatrixXd numericalJacobian(
        const Eigen::VectorXd& extendedState, const Eigen::VectorXd& Hx, const Eigen::VectorXd& F0
    ) const;

    /**
     * @brief Compute the sign-oriented unit tangent to ker(A) (A is n x (n+1), rank n), via QR of
     *        A^T, along with the QR factors used by the corrector step. Port of numcontin.m's
     *        static tangent().
     */
    struct TangentResult
    {
        Eigen::VectorXd tangent;
        Eigen::MatrixXd Q;
        Eigen::MatrixXd R;
    };
    static TangentResult tangent(const Eigen::MatrixXd& A, double orient);
};
