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

#include "../src/numerics/MCSCContinuationSolver.h"

#include <gtest/gtest.h>
#include <cmath>

// Synthetic root-finding problems, independent of any MCSC domain machinery, so continuation
// algorithm bugs are distinguishable from MCSC-residual bugs.

TEST(MCSCContinuationSolverTest, SolvesLinearSystem)
{
    // F(X) = M*X - b, a trivial linear system with a known unique root.
    Eigen::MatrixXd M(2, 2);
    M << 2.0, 1.0, 1.0, 3.0;
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    auto F = [&](const Eigen::VectorXd& X) -> Eigen::VectorXd { return M * X - b; };

    Eigen::VectorXd x0(2);
    x0 << 0.0, 0.0;

    MCSCContinuationSolver solver(F);
    auto result = solver.solve(x0);

    Eigen::VectorXd expected = M.colPivHouseholderQr().solve(b);
    EXPECT_NEAR(result.solution[0], expected[0], 1e-8);
    EXPECT_NEAR(result.solution[1], expected[1], 1e-8);
    EXPECT_LT(result.residualNorm, 1e-10);
    EXPECT_NEAR(result.lastLambda, 1.0, 1e-9);
}

TEST(MCSCContinuationSolverTest, SolvesMildlyNonlinearSystem)
{
    // F1 = x^2 + y^2 - 5, F2 = x - y - 1: root at (x,y) = (2,1) (also (-1,-2), but continuation
    // from a positive-quadrant initial guess should track to (2,1)).
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd
    {
        Eigen::VectorXd result(2);
        result[0] = X[0] * X[0] + X[1] * X[1] - 5.0;
        result[1] = X[0] - X[1] - 1.0;
        return result;
    };

    Eigen::VectorXd x0(2);
    x0 << 3.0, 3.0;

    MCSCContinuationSolver solver(F);
    auto result = solver.solve(x0);

    EXPECT_NEAR(result.solution[0], 2.0, 1e-6);
    EXPECT_NEAR(result.solution[1], 1.0, 1e-6);
    EXPECT_LT(result.residualNorm, 1e-10);
}

TEST(MCSCContinuationSolverTest, ResultIsAlreadyASolutionWhenInitialGuessIsExact)
{
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd
    {
        Eigen::VectorXd result(1);
        result[0] = X[0] - 3.0;
        return result;
    };

    Eigen::VectorXd x0(1);
    x0 << 3.0;

    MCSCContinuationSolver solver(F);
    auto result = solver.solve(x0);

    EXPECT_NEAR(result.solution[0], 3.0, 1e-9);
    EXPECT_LT(result.residualNorm, 1e-10);
}

TEST(MCSCContinuationSolverTest, ThrowsWhenMaxStepsExceeded)
{
    // A function whose homotopy never settles below tolerance within a tiny step budget --
    // forces the maxSteps safety cap rather than the hmin path.
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd
    {
        Eigen::VectorXd result(1);
        result[0] = std::sin(1000.0 * X[0]);
        return result;
    };

    Eigen::VectorXd x0(1);
    x0 << 0.6;

    MCSCContinuationSolver::Options opts;
    opts.maxSteps = 2;
    MCSCContinuationSolver solver(F, opts);

    try
    {
        solver.solve(x0);
        FAIL() << "expected ConvergenceError";
    }
    catch (const MCSCContinuationSolver::ConvergenceError& e)
    {
        EXPECT_NE(std::string(e.what()).find("maximum number of continuation steps"), std::string::npos);
    }
}

TEST(MCSCContinuationSolverTest, ThrowsOnHminUnderflowBeforeReachingTolerance)
{
    // An oscillatory homotopy (like ThrowsWhenMaxStepsExceeded's) repeatedly fails the
    // angle/contraction checks and halves h every time, but with maxSteps left generous the step
    // count never triggers first -- h itself underflows hmin, exercising a distinct throw path.
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd
    {
        Eigen::VectorXd result(1);
        result[0] = std::sin(1000.0 * X[0]);
        return result;
    };

    Eigen::VectorXd x0(1);
    x0 << 0.6;

    MCSCContinuationSolver::Options opts;
    opts.maxSteps = 100000;
    MCSCContinuationSolver solver(F, opts);

    try
    {
        solver.solve(x0);
        FAIL() << "expected ConvergenceError";
    }
    catch (const MCSCContinuationSolver::ConvergenceError& e)
    {
        EXPECT_NE(std::string(e.what()).find("step size underflowed hmin"), std::string::npos);
    }
}

TEST(MCSCContinuationSolverTest, ThrowsOnNonFiniteResidual)
{
    // F blows up away from x=0, so continuation stepping past that point produces a non-finite
    // residual -- verifies the explicit isfinite guard rather than relying on IEEE 754 comparison
    // semantics (NaN > tolerance is false) to silently exit as "converged".
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd
    {
        Eigen::VectorXd result(1);
        result[0] = 1.0 / X[0];
        return result;
    };

    Eigen::VectorXd x0(1);
    x0 << 1.0;

    MCSCContinuationSolver::Options opts;
    opts.hmax = 100.0;
    opts.initialStep = 50.0;
    MCSCContinuationSolver solver(F, opts);

    EXPECT_THROW(solver.solve(x0), MCSCContinuationSolver::ConvergenceError);
}

TEST(MCSCContinuationSolverTest, ConstructorThrowsOnInvalidOptions)
{
    auto F = [](const Eigen::VectorXd& X) -> Eigen::VectorXd { return X; };

    MCSCContinuationSolver::Options negativeTolerance;
    negativeTolerance.tolerance = -1e-10;
    EXPECT_THROW(MCSCContinuationSolver(F, negativeTolerance), std::invalid_argument);

    MCSCContinuationSolver::Options zeroHmin;
    zeroHmin.hmin = 0.0;
    EXPECT_THROW(MCSCContinuationSolver(F, zeroHmin), std::invalid_argument);

    MCSCContinuationSolver::Options hmaxNotGreaterThanHmin;
    hmaxNotGreaterThanHmin.hmax = hmaxNotGreaterThanHmin.hmin;
    EXPECT_THROW(MCSCContinuationSolver(F, hmaxNotGreaterThanHmin), std::invalid_argument);

    MCSCContinuationSolver::Options zeroFiniteDifferenceStep;
    zeroFiniteDifferenceStep.finiteDifferenceStep = 0.0;
    EXPECT_THROW(MCSCContinuationSolver(F, zeroFiniteDifferenceStep), std::invalid_argument);

    MCSCContinuationSolver::Options zeroMaxSteps;
    zeroMaxSteps.maxSteps = 0;
    EXPECT_THROW(MCSCContinuationSolver(F, zeroMaxSteps), std::invalid_argument);
}
