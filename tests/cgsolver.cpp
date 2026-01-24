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