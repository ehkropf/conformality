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

#include <gtest/gtest.h>
#include "../src/numerics/CGSolver.h"
#include "../src/methods/FornbergMCConfiguration.h"
#include "../src/core/StatusManager.h"
#include "../src/core/Types.h"

class CGSolverTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create test configuration
        config.N = 64; // Smaller N for faster tests
        config.max_newton_iterations = 10;
        config.max_cgm_iterations = 100;
        config.newton_tolerance = 1e-8;
        config.cgm_tolerance = 1e-8;
        config.verbose = false; // Quiet tests
    }

    FornbergMCConfiguration config;
};

// Test CGSolver
TEST_F(CGSolverTest, Construction)
{
    // Construction should work
    EXPECT_NO_THROW(CGSolver solver(config));
}

TEST_F(CGSolverTest, SelfTest)
{
    CGSolver solver(config);

    // Self-test should pass
    EXPECT_TRUE(solver.runSelfTest(50));
}

TEST_F(CGSolverTest, SimpleSystemSolution)
{
    CGSolver solver(config);
    
    // Test simple system solution
    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    
    // Identity matrix function: A*x = x
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };
    
    Eigen::VectorXd solution = solver.solve(identity_function, b);
    
    // Solution should be close to b for identity system
    EXPECT_LT((solution - b).norm(), 1e-10);
    EXPECT_TRUE(solver.hasConverged());
    EXPECT_LT(solver.getFinalResidual(), config.cgm_tolerance);
}

TEST_F(CGSolverTest, ConvergenceInfo)
{
    CGSolver solver(config);
    
    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };
    
    solver.solve(identity_function, b);
    
    // Check convergence information
    const auto& info = solver.getLastConvergenceInfo();
    EXPECT_TRUE(info.converged);
    EXPECT_GT(info.iterations, 0);
    EXPECT_LT(info.final_residual, config.cgm_tolerance);
    EXPECT_FALSE(info.residual_history.empty());
}

TEST_F(CGSolverTest, ComplexSystemSolution)
{
    CGSolver solver(config);
    
    int n = 8;
    
    // Simple complex system: D*x = [1+i, 2+2i, ...]
    auto D_function = [n](const Eigen::VectorXd& x) -> Eigen::VectorXcd {
        Eigen::VectorXcd result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = Complex(x[i], x[i]); // (1+i) * x[i]
        }
        return result;
    };
    
    auto D_adjoint_function = [n](const Eigen::VectorXcd& x) -> Eigen::VectorXcd {
        Eigen::VectorXcd result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = Complex(x[i].real() - x[i].imag(), x[i].real() + x[i].imag()); // (1-i) * x[i]
        }
        return result;
    };
    
    Eigen::VectorXcd g = Eigen::VectorXcd::Ones(n);
    
    // This should solve the system without throwing
    EXPECT_NO_THROW(solver.solveComplexSystem(D_function, D_adjoint_function, g));
}

TEST_F(CGSolverTest, ComplexSystemConvergenceInfo)
{
    CGSolver solver(config);
    
    int n = 6;
    
    // Simple diagonal complex system
    auto D_function = [n](const Eigen::VectorXd& x) -> Eigen::VectorXcd {
        Eigen::VectorXcd result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = Complex(2.0 * x[i], x[i]); // (2+i) * x[i]
        }
        return result;
    };
    
    auto D_adjoint_function = [n](const Eigen::VectorXcd& x) -> Eigen::VectorXcd {
        Eigen::VectorXcd result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = Complex(2.0 * x[i].real() - x[i].imag(), 2.0 * x[i].imag() + x[i].real()); // (2-i) * x[i]
        }
        return result;
    };
    
    Eigen::VectorXcd g = Eigen::VectorXcd::Ones(n);
    
    Eigen::VectorXd solution = solver.solveComplexSystem(D_function, D_adjoint_function, g);
    
    // Check that we got a reasonable solution
    EXPECT_EQ(solution.size(), n);
    EXPECT_TRUE(solver.hasConverged());
    
    const auto& info = solver.getLastConvergenceInfo();
    EXPECT_TRUE(info.converged);
    EXPECT_GT(info.iterations, 0);
}

TEST_F(CGSolverTest, StatusManagerLogging)
{
    auto statusManager = std::make_shared<StatusManager>();
    CGSolver solver(config);
    solver.setStatusManager(statusManager);

    // Run a simple solve which should generate log messages
    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    solver.solve(identity_function, b);

    // Check that log messages were generated
    auto messages = statusManager->getMessages();
    EXPECT_FALSE(messages.empty()) << "Expected log messages from CGSolver";

    // Should have at least debug messages for starting CG and convergence info
    auto debugMessages = statusManager->getMessages(StatusLevel::DEBUG);
    EXPECT_FALSE(debugMessages.empty()) << "Expected DEBUG messages from CGSolver";

    // Check for expected component name in messages
    bool found_cgsolver_message = false;
    for (const auto& msg : messages)
    {
        if (msg.component == "CGSolver")
        {
            found_cgsolver_message = true;
            break;
        }
    }
    EXPECT_TRUE(found_cgsolver_message) << "Expected messages with component 'CGSolver'";

    // Should have an INFO message for successful convergence
    auto infoMessages = statusManager->getMessages(StatusLevel::INFO);
    EXPECT_FALSE(infoMessages.empty()) << "Expected INFO message for convergence";
}

TEST_F(CGSolverTest, StatusManagerSelfTestLogging)
{
    auto statusManager = std::make_shared<StatusManager>();
    CGSolver solver(config);
    solver.setStatusManager(statusManager);

    // Run self-test which should log messages
    bool result = solver.runSelfTest(50);
    EXPECT_TRUE(result);

    // Check that log messages were generated
    auto messages = statusManager->getMessages();
    EXPECT_FALSE(messages.empty()) << "Expected log messages from CGSolver self-test";

    // Check for debug message about starting self-test
    bool found_selftest_start = false;
    for (const auto& msg : messages)
    {
        if (msg.message.find("self-test") != std::string::npos ||
            msg.message.find("Self-test") != std::string::npos)
        {
            found_selftest_start = true;
            break;
        }
    }
    EXPECT_TRUE(found_selftest_start) << "Expected self-test message";
}

TEST_F(CGSolverTest, NoLoggingWithoutStatusManager)
{
    // Without a StatusManager, solve should still work
    CGSolver solver(config);
    // Do NOT set status manager

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    // Should not throw and should still converge
    EXPECT_NO_THROW(solver.solve(identity_function, b));
    EXPECT_TRUE(solver.hasConverged());
}

// GH-92: CGSolver should exit early on non-finite input instead of looping

TEST_F(CGSolverTest, EarlyExitOnInfRHS)
{
    CGSolver solver(config);
    solver.setStatusManager(std::make_shared<StatusManager>());

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    b(0) = std::numeric_limits<double>::infinity();
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    solver.solve(identity_function, b);
    EXPECT_FALSE(solver.hasConverged());
    EXPECT_EQ(solver.getLastConvergenceInfo().iterations, 0);
}

TEST_F(CGSolverTest, EarlyExitOnNanRHS)
{
    CGSolver solver(config);
    solver.setStatusManager(std::make_shared<StatusManager>());

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    b(3) = std::numeric_limits<double>::quiet_NaN();
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    solver.solve(identity_function, b);
    EXPECT_FALSE(solver.hasConverged());
    EXPECT_EQ(solver.getLastConvergenceInfo().iterations, 0);
}

TEST_F(CGSolverTest, EarlyExitOnNonFiniteMatVecProduct)
{
    CGSolver solver(config);
    solver.setStatusManager(std::make_shared<StatusManager>());

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    // Matrix-vector product that always returns a vector with inf in position 0
    auto bad_function = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd result = x;
        result(0) = std::numeric_limits<double>::infinity();
        return result;
    };

    solver.solve(bad_function, b);
    // Should not converge, and should not run the full max_cgm_iterations
    EXPECT_FALSE(solver.hasConverged());
    EXPECT_LT(solver.getLastConvergenceInfo().iterations, config.max_cgm_iterations);
}

TEST_F(CGSolverTest, ThrowsOnInfRHSWithoutStatusManager)
{
    CGSolver solver(config);
    // Do NOT set StatusManager

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    b(0) = std::numeric_limits<double>::infinity();
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    EXPECT_THROW(solver.solve(identity_function, b), std::runtime_error);
}

TEST_F(CGSolverTest, ThrowsOnNanRHSWithoutStatusManager)
{
    CGSolver solver(config);
    // Do NOT set StatusManager

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    b(3) = std::numeric_limits<double>::quiet_NaN();
    auto identity_function = [](const Eigen::VectorXd& x) { return x; };

    EXPECT_THROW(solver.solve(identity_function, b), std::runtime_error);
}

TEST_F(CGSolverTest, ThrowsOnNonFiniteMatVecProductWithoutStatusManager)
{
    CGSolver solver(config);
    // Do NOT set StatusManager

    int n = 10;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    auto bad_function = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd result = x;
        result(0) = std::numeric_limits<double>::infinity();
        return result;
    };

    EXPECT_THROW(solver.solve(bad_function, b), std::runtime_error);
}
