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

#include "FornbergMCConfiguration.h"
#include <stdexcept>

void FornbergMCConfiguration::validate() const
{
    // NOTE: This validation covers parameter ranges and logical consistency,
    // but does NOT validate:
    // - Algorithm-performance relationships (e.g., max_cgm_iterations scaling with N)
    // - Domain-dependent parameter appropriateness (e.g., N sufficient for domain complexity)
    // - Hardware constraints (e.g., max_N fitting in available memory)
    // - Runtime context dependencies (e.g., minimum_progress_rate for specific problems)
    // - Cross-parameter optimization (e.g., tolerance/iteration limit pairing)

    // Newton iteration parameters
    if (newton_tolerance <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: newton_tolerance must be positive");
    }
    if (newton_tolerance < 1e-15)
    {
        throw std::invalid_argument("FornbergMCConfiguration: newton_tolerance too small, approaching machine precision");
    }

    if (max_newton_iterations <= 0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: max_newton_iterations must be positive");
    }

    if (newton_damping_factor <= 0.0 || newton_damping_factor > 1.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: newton_damping_factor must be in (0, 1]");
    }

    // Conjugate gradient parameters
    if (cgm_tolerance <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: cgm_tolerance must be positive");
    }
    if (cgm_tolerance < 1e-15)
    {
        throw std::invalid_argument("FornbergMCConfiguration: cgm_tolerance too small, approaching machine precision");
    }

    if (max_cgm_iterations <= 0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: max_cgm_iterations must be positive");
    }

    if (cgm_restart_threshold <= 0.0 || cgm_restart_threshold >= 1.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: cgm_restart_threshold must be in (0, 1)");
    }

    if (max_cgm_restarts < 0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: max_cgm_restarts must be non-negative");
    }

    // Tolerance consistency check
    if (cgm_tolerance > newton_tolerance)
    {
        throw std::invalid_argument("FornbergMCConfiguration: cgm_tolerance should be <= newton_tolerance for efficiency");
    }

    // Matrix conditioning parameters
    if (eigenvalue_warning_threshold <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: eigenvalue_warning_threshold must be positive");
    }

    if (enable_regularization && regularization_parameter <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: regularization_parameter must be positive when enabled");
    }

    // Discretization parameters
    if (N <= 0 || (N & (N - 1)) != 0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: N must be a positive power of 2");
    }

    if (max_N < N)
    {
        throw std::invalid_argument("FornbergMCConfiguration: max_N must be >= N");
    }

    if (adaptive_N && (max_N <= 0 || (max_N & (max_N - 1)) != 0))
    {
        throw std::invalid_argument("FornbergMCConfiguration: max_N must be a positive power of 2 when adaptive_N is enabled");
    }

    if (series_truncation_tolerance <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: series_truncation_tolerance must be positive");
    }

    // Redistribution parameters
    if (redistribution_threshold <= 0.0 || redistribution_threshold > 1.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: redistribution_threshold must be in (0, 1]");
    }

    if (redistribution_frequency <= 0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: redistribution_frequency must be positive");
    }

    // Annulus parameters
    if (annulus_aspect_ratio_threshold <= 1.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: annulus_aspect_ratio_threshold must be > 1");
    }

    // Error handling parameters
    if (minimum_progress_rate <= 0.0)
    {
        throw std::invalid_argument("FornbergMCConfiguration: minimum_progress_rate must be positive");
    }

    // Performance parameters
    if (parallel_threads < 1)
    {
        throw std::invalid_argument("FornbergMCConfiguration: parallel_threads must be >= 1");
    }
}

void FornbergMCConfiguration::setAnnulusOptimized()
{
    auto_detect_annulus = true;
    force_annulus_mode = false;
    newton_tolerance = 1e-14;
    cgm_tolerance = 1e-14;
    enable_best_iterate = true;
}

void FornbergMCConfiguration::setHighPrecision()
{
    newton_tolerance = 1e-15;
    cgm_tolerance = 1e-15;
    series_truncation_tolerance = 1e-16;
    enable_regularization = false;
    monitor_eigenvalues = true;
}

void FornbergMCConfiguration::setFastComputation()
{
    newton_tolerance = 1e-8;
    cgm_tolerance = 1e-8;
    max_newton_iterations = 20;
    max_cgm_iterations = 20;
    max_cgm_restarts = 0;
    verbose = false;
    boundary_quality_check = false;
}
