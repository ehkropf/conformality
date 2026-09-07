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

#include "MCSCContinuationSolver.h"

#include "../core/Tolerances.h"

#include <algorithm>
#include <cmath>

void MCSCContinuationSolver::validateOptions(const Options& options)
{
    if (options.tolerance <= 0.0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: tolerance must be positive");
    }
    if (options.hmin <= 0.0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: hmin must be positive");
    }
    if (options.hmax <= options.hmin)
    {
        throw std::invalid_argument("MCSCContinuationSolver: hmax must be greater than hmin");
    }
    if (options.finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: finiteDifferenceStep must be positive");
    }
    if (options.initialStep <= 0.0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: initialStep must be positive");
    }
    if (options.maxSteps <= 0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: maxSteps must be positive");
    }
    if (options.maxRefinementsBeforeJacobianReset < 0)
    {
        throw std::invalid_argument("MCSCContinuationSolver: maxRefinementsBeforeJacobianReset must be non-negative");
    }
}

MCSCContinuationSolver::MCSCContinuationSolver(ResidualFunction F)
    : MCSCContinuationSolver(std::move(F), Options{})
{
}

MCSCContinuationSolver::MCSCContinuationSolver(ResidualFunction F, Options options)
    : m_F{std::move(F)}
    , m_options{options}
{
    validateOptions(m_options);
}

Eigen::VectorXd MCSCContinuationSolver::curveFunction(
    const Eigen::VectorXd& extendedState, const Eigen::VectorXd& F0
) const
{
    const double lambda = extendedState[0];
    const Eigen::VectorXd X = extendedState.tail(extendedState.size() - 1);
    return m_F(X) + (lambda - 1.0) * F0;
}

Eigen::MatrixXd MCSCContinuationSolver::numericalJacobian(
    const Eigen::VectorXd& extendedState, const Eigen::VectorXd& Hx, const Eigen::VectorXd& F0
) const
{
    const int n = static_cast<int>(Hx.size());
    Eigen::MatrixXd J(n, n + 1);
    const double h = m_options.finiteDifferenceStep;

    for (int k = 0; k < n + 1; ++k)
    {
        Eigen::VectorXd perturbed = extendedState;
        perturbed[k] += h;
        J.col(k) = (curveFunction(perturbed, F0) - Hx) / h;
    }

    return J;
}

MCSCContinuationSolver::TangentResult MCSCContinuationSolver::tangent(const Eigen::MatrixXd& A, double orient)
{
    // A is n x (n+1); qr(A') in MATLAB is the full QR of A^T ((n+1) x n), giving an (n+1) x (n+1)
    // orthogonal Q whose last column spans ker(A) (rank(A) == n assumed, matching numcontin.m's
    // implicit assumption of a smooth curve with nowhere-vanishing tangent).
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A.transpose());
    const int n = static_cast<int>(A.rows());
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n + 1, n + 1);
    Eigen::MatrixXd R = qr.matrixQR().topRows(n).template triangularView<Eigen::Upper>();

    // A near-zero diagonal entry means A has become (nearly) rank-deficient -- the curve's
    // tangent direction is no longer well-defined, and the sign computation below would be
    // dictated by numerical noise rather than genuine orientation. This is exactly the kind of
    // silently-degrading condition CLAUDE.md calls out as needing to surface, not be absorbed
    // into a plausible-looking but meaningless step.
    const double maxDiag = R.diagonal().cwiseAbs().maxCoeff();
    const double minDiag = R.diagonal().cwiseAbs().minCoeff();
    if (maxDiag > 0.0 && minDiag < 1e-10 * maxDiag)
    {
        throw ConvergenceError("MCSCContinuationSolver: continuation Jacobian became rank-deficient");
    }

    Eigen::VectorXd t = Q.col(n);

    double detQ = Q.determinant();
    double signProd = 1.0;
    for (int i = 0; i < n; ++i)
    {
        signProd *= (R(i, i) >= 0.0 ? 1.0 : -1.0);
    }
    const double sign = orient * (detQ >= 0.0 ? 1.0 : -1.0) * signProd;

    TangentResult result;
    result.tangent = sign * t;
    result.Q = std::move(Q);
    result.R = std::move(R);
    return result;
}

MCSCContinuationSolver::Result MCSCContinuationSolver::solve(const Eigen::VectorXd& x0) const
{
    const int n = static_cast<int>(x0.size());
    const Eigen::VectorXd F0 = m_F(x0);

    Eigen::VectorXd x(n + 1);
    x[0] = 0.0;
    x.tail(n) = x0;

    Eigen::VectorXd refdir = Eigen::VectorXd::Zero(n + 1);
    refdir[0] = 1.0;

    Eigen::VectorXd Fx = curveFunction(x, F0);
    Eigen::MatrixXd A = numericalJacobian(x, Fx, F0);

    TangentResult tr = tangent(A, 1.0);
    Eigen::VectorXd t = tr.tangent;
    const double orient = (t.dot(refdir) >= 0.0) ? 1.0 : -1.0;
    t *= orient;

    double h = m_options.initialStep;
    bool adaptNewton = false;
    int refineCount = 0;
    double objectiveNorm = F0.lpNorm<Eigen::Infinity>();
    // Seed the loop guard so an exact initial root (F0 == 0, i.e. x0 is already a solution of the
    // original F) does not read as "already converged" before a single continuation step has
    // actually verified lambda reaches 1 -- overwritten by the true homotopy residual on the
    // first iteration below.
    if (objectiveNorm == 0.0)
    {
        objectiveNorm = 1.0;
    }

    int steps = 0;
    while (std::abs(h) > m_options.hmin && objectiveNorm > m_options.tolerance)
    {
        if (steps++ > m_options.maxSteps)
        {
            throw ConvergenceError(
                "MCSCContinuationSolver: exceeded maximum number of continuation steps");
        }

        if (refineCount > m_options.maxRefinementsBeforeJacobianReset)
        {
            A = numericalJacobian(x, Fx, F0);
            tr = tangent(A, orient);
            t = tr.tangent;
            refineCount = 0;
        }

        const Eigen::VectorXd u = x + h * t;
        const Eigen::VectorXd Fu = curveFunction(u, F0);

        // Predictor update: Broyden rank-one update to the Jacobian, then re-derive the tangent.
        const Eigen::VectorXd s = t;
        A = A + (Fu - Fx) * t.transpose() / h;
        tr = tangent(A, orient);
        t = tr.tangent;

        const double cosAngle = std::clamp(s.dot(t), -1.0, 1.0);
        if (std::acos(cosAngle) > m_options.angleMax)
        {
            h /= 2.0;
            ++refineCount;
            continue;
        }

        // Corrector step: solve via the QR factors from the (updated) tangent's Jacobian.
        const Eigen::VectorXd du = tr.Q.leftCols(n)
            * (tr.R.topLeftCorner(n, n).transpose().triangularView<Eigen::Lower>().solve(Fu));
        const Eigen::VectorXd v = u - du;
        const Eigen::VectorXd Fv = curveFunction(v, F0);

        // Corrector update: Broyden rank-one update using the corrector displacement. A
        // denormal-but-nonzero du would otherwise blow up this update into a meaningless
        // Jacobian correction, so the guard uses PIVOT_EPS (a divisor-magnitude floor), not a
        // bare nonzero check.
        const double duNormSq = du.squaredNorm();
        if (duNormSq > PIVOT_EPS)
        {
            A = A + Fv * (-du).transpose() / duNormSq;
        }
        tr = tangent(A, orient);
        t = tr.tangent;

        const double contract = Fv.norm() / (Fu.norm() + m_options.resmin);
        if (contract > m_options.contractionMax)
        {
            h /= 2.0;
            ++refineCount;
            continue;
        }

        objectiveNorm = (Fv - (v[0] - 1.0) * F0).lpNorm<Eigen::Infinity>();
        // A non-finite residual (e.g. F evaluated a degenerate configuration mid-step) must be
        // caught explicitly here: under IEEE 754 semantics, any comparison against NaN is false,
        // so both the while-loop guard above and the post-loop tolerance check below would
        // otherwise treat a NaN residual as "already converged" and silently return a meaningless
        // solution instead of throwing.
        if (!std::isfinite(objectiveNorm))
        {
            throw ConvergenceError(
                "MCSCContinuationSolver: encountered non-finite residual during continuation (step "
                + std::to_string(steps) + ")");
        }

        x = v;
        Fx = Fv;
        refineCount = 0;

        constexpr double target = 1.0;
        if (v[0] > target)
        {
            adaptNewton = true;
        }

        if (adaptNewton)
        {
            if (std::abs(t[0]) < m_options.hmin)
            {
                throw ConvergenceError(
                    "MCSCContinuationSolver: tangent's lambda-component vanished during final Newton phase");
            }
            h = -(v[0] - target) / t[0];
        }
        else
        {
            h = std::min(m_options.hmax, 2.0 * h);
        }
    }

    if (objectiveNorm > m_options.tolerance)
    {
        throw ConvergenceError(
            "MCSCContinuationSolver: step size underflowed hmin before reaching tolerance "
            "(residual norm " + std::to_string(objectiveNorm) + ")");
    }

    Result result;
    result.solution = x.tail(n);
    result.residualNorm = objectiveNorm;
    result.lastLambda = x[0];
    result.lastStepSize = h;
    return result;
}
