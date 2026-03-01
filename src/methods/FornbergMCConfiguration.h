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

/**
 * @brief Configuration structure for the Fornberg-like method for multiply connected domains
 *
 * This structure contains all parameters and settings for the FornbergMC implementation,
 * including Newton iteration parameters, conjugate gradient solver settings, matrix 
 * conditioning parameters, boundary parameterization control, and diagnostic options.
 */
struct FornbergMCConfiguration
{
    // Newton iteration parameters
    double newton_tolerance = 1e-12;        // Convergence threshold for ||update||
    int max_newton_iterations = 50;         // Maximum Newton steps before failure
    bool enable_newton_damping = false;     // MATLAB uses undamped Newton (Phase 1 parity). Consider enabling in Phase 2 for robustness.
    double newton_damping_factor = 0.5;     // Step reduction multiplier when damping enabled
    
    // Conjugate gradient solver parameters  
    double cgm_tolerance = 1e-12;           // Relative residual tolerance for CG
    int max_cgm_iterations = 1000;          // Maximum CG iterations (should scale with N)
    bool enable_best_iterate = true;        // Track best iterate for non-convergent cases
    double cgm_restart_threshold = 0.1;     // Restart CG if residual stagnates
    int max_cgm_restarts = 5;               // Maximum CG restarts before giving up (0 disables restarts)
    
    // Matrix conditioning and stability
    bool monitor_eigenvalues = true;        // Check for nullspace issues
    double eigenvalue_warning_threshold = 1e-10;  // Warn if eigenvalues approach machine epsilon
    bool enable_regularization = false;     // Add small diagonal term if ill-conditioned
    double regularization_parameter = 1e-14; // Regularization strength
    
    // Boundary parameterization control
    bool enable_redistribution = true;      // Automatic boundary point redistribution
    double redistribution_threshold = 0.1;  // Trigger when points cluster/separate
    int redistribution_frequency = 5;       // Check every N Newton iterations
    bool preserve_boundary_orientation = true; // Maintain counterclockwise ordering
    
    // Numerical discretization parameters
    int N = 128;                           // Boundary points per component (power of 2)
    bool adaptive_N = false;               // Automatically increase N if needed
    int max_N = 512;                       // Upper limit for adaptive refinement
    double series_truncation_tolerance = 1e-14; // Drop high-frequency coefficients
    
    // Initial guess strategy
    enum class InitialGuessMethod
    {
        MANUAL,                            // User-provided c_ν, ρ_ν values
        GEOMETRIC_HEURISTIC,               // Centroid + inscribed circle radius
        BOUNDARY_ANALYSIS,                 // Analyze curvature and convex hull
        PREVIOUS_SOLUTION                  // Continuation from similar domain
    };
    InitialGuessMethod initial_guess_method = InitialGuessMethod::GEOMETRIC_HEURISTIC;
    
    // Annulus case optimization (m=2)
    bool auto_detect_annulus = true;       // Automatically use specialized formulation
    bool force_annulus_mode = false;       // Override detection for testing
    double annulus_aspect_ratio_threshold = 10.0; // Switch to general method if extreme
    
    // Diagnostics and monitoring
    bool verbose = false;                  // Print iteration progress
    bool save_iteration_history = false;   // Store Newton updates for analysis
    bool eigenvalue_analysis = false;      // Compute and display matrix spectrum
    bool boundary_quality_check = true;    // Validate parameterization quality
    
    // Error handling and recovery
    bool enable_fallback_methods = true;   // Try alternative approaches on failure
    bool strict_convergence = false;       // Fail if any tolerance not met
    double minimum_progress_rate = 1e-6;   // Minimum Newton update norm reduction
    
    // Performance tuning
    bool use_fft_optimization = true;      // Optimize FFT calls for repeated evaluation
    bool cache_matrix_operations = false;  // Cache expensive matrix computations
    int parallel_threads = 1;             // Thread count for parallel operations (future)

    /**
     * @brief Validate the configuration parameters
     * @throws std::invalid_argument if any parameter is invalid
     */
    void validate() const;

    /**
     * @brief Set parameters for annulus-optimized computation
     */
    void setAnnulusOptimized();

    /**
     * @brief Set parameters for high-precision computation
     */
    void setHighPrecision();

    /**
     * @brief Set parameters for fast computation (lower precision)
     */
    void setFastComputation();
};
