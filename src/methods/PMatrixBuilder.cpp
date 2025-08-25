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

#include "PMatrixBuilder.h"
// #include "../core/StatusManager.h"
#include <stdexcept>

PMatrixBuilder::PMatrixBuilder(const FornbergMCConfiguration& config, int connectivity, bool is_annulus)
    : m_config{config}
    , m_connectivity{connectivity}
    , m_N{config.N}
    , m_is_annulus{is_annulus}
{
    validateParameters();
    
    // Initialize frequency indices and normalization conditions
    initializeFrequencyIndices();
    setupNormalizationConditions();
    
    // Logging removed for compilation
}

Eigen::MatrixXcd PMatrixBuilder::buildPMatrix(int nu) const
{
    if (nu < 0 || nu >= m_connectivity)
    {
        throw std::invalid_argument("PMatrixBuilder: Invalid component index " + std::to_string(nu));
    }
    
    if (m_is_annulus)
    {
        return buildAnnulusPMatrix(nu);
    }
    else
    {
        return buildGeneralPMatrix(nu);
    }
}

std::vector<Eigen::MatrixXcd> PMatrixBuilder::buildAllPMatrices() const
{
    std::vector<Eigen::MatrixXcd> P_matrices;
    P_matrices.reserve(m_connectivity);
    
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        P_matrices.push_back(buildPMatrix(nu));
    }
    
    // Debug logging removed for compilation
    return P_matrices;
}

const Eigen::VectorXi& PMatrixBuilder::getFrequencyIndices(int nu) const
{
    if (nu < 0 || nu >= m_connectivity)
    {
        throw std::invalid_argument("PMatrixBuilder: Invalid component index " + std::to_string(nu));
    }
    
    return m_frequency_indices[nu];
}

void PMatrixBuilder::applyNormalizationConditions(Eigen::MatrixXcd& system_matrix, 
                                                Eigen::VectorXcd& rhs_vector) const
{
    // Debug logging removed for compilation
    
    // TODO: Implement normalization condition application
    // This involves modifying specific rows/columns of the system matrix
    // to enforce the 3 degrees of freedom normalization per Newton step
    
    for (const auto& condition : m_normalization_conditions)
    {
        // Placeholder: modify system matrix to enforce normalization
        // The actual implementation depends on the specific normalization strategy
        // Debug logging removed for compilation
    }
}

int PMatrixBuilder::getSystemSize() const
{
    // Total system size depends on connectivity and discretization
    // For m boundary components with N points each, the system size is approximately:
    // - General case: sum of analyticity constraints for each component
    // - Annulus case: optimized size due to specialized structure
    
    if (m_is_annulus)
    {
        // Annulus case has specialized structure
        return 2 * m_N - 6; // Approximate size for annulus
    }
    else
    {
        // General case: each component contributes constraints
        return m_connectivity * m_N - 3 * (m_connectivity - 1); // 3m-6 moduli parameters
    }
}

void PMatrixBuilder::setConnectivity(int new_connectivity)
{
    if (new_connectivity < 2)
    {
        throw std::invalid_argument("PMatrixBuilder: Connectivity must be at least 2");
    }
    
    m_connectivity = new_connectivity;
    m_is_annulus = (new_connectivity == 2) && m_config.auto_detect_annulus;
    
    // Reinitialize with new connectivity
    initializeFrequencyIndices();
    setupNormalizationConditions();
    
    // Logging removed for compilation
}

void PMatrixBuilder::setAnnulusMode(bool use_annulus)
{
    if (use_annulus && m_connectivity != 2)
    {
        // Warning logging removed for compilation
        return;
    }
    
    m_is_annulus = use_annulus;
    
    // Reinitialize frequency indices for new mode
    initializeFrequencyIndices();
    setupNormalizationConditions();
    
    // Logging removed for compilation
}

void PMatrixBuilder::initializeFrequencyIndices()
{
    m_frequency_indices.clear();
    m_frequency_indices.resize(m_connectivity);
    
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        // Initialize frequency indices based on FFT structure
        // FFTW output order: [0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]
        
        m_frequency_indices[nu].resize(m_N);
        for (int j = 0; j < m_N; ++j)
        {
            if (j <= m_N/2)
            {
                m_frequency_indices[nu][j] = j;  // Positive frequencies and zero
            }
            else
            {
                m_frequency_indices[nu][j] = j - m_N;  // Negative frequencies
            }
        }
    }
    
    // Debug logging removed for compilation
}

Eigen::MatrixXcd PMatrixBuilder::buildGeneralPMatrix(int nu) const
{
    // Stub implementation for general P matrix construction
    // Debug logging removed for compilation
    
    // Determine matrix size based on analyticity conditions
    int positive_freq_count = getPositiveFrequencyCount(nu);
    int total_size = m_N;
    
    Eigen::MatrixXcd P_nu = Eigen::MatrixXcd::Zero(positive_freq_count, total_size);
    
    // TODO: Implement general P matrix construction
    // This involves:
    // 1. Identifying positive frequency coefficients to zero
    // 2. Setting up linear constraints to enforce analyticity
    // 3. Handling boundary component-specific frequency structure
    
    // Placeholder: create identity-like structure for positive frequencies
    for (int i = 0; i < positive_freq_count && i < total_size; ++i)
    {
        P_nu(i, i + 1) = 1.0;  // Zero positive frequencies 1, 2, ..., N/2-1
    }
    
    return P_nu;
}

Eigen::MatrixXcd PMatrixBuilder::buildAnnulusPMatrix(int nu) const
{
    // Stub implementation for annulus-optimized P matrix construction
    // Debug logging removed for compilation
    
    // Annulus case has specialized structure for m=2
    if (m_connectivity != 2)
    {
        throw std::runtime_error("PMatrixBuilder: Annulus P matrix requested but connectivity is not 2");
    }
    
    int positive_freq_count = getPositiveFrequencyCount(nu);
    int total_size = m_N;
    
    Eigen::MatrixXcd P_nu = Eigen::MatrixXcd::Zero(positive_freq_count, total_size);
    
    // TODO: Implement annulus-optimized P matrix construction
    // This involves specialized handling for the 2-connected case
    // with optimized constraint structure
    
    // Placeholder: similar to general case but optimized
    for (int i = 0; i < positive_freq_count && i < total_size; ++i)
    {
        P_nu(i, i + 1) = 1.0;
    }
    
    return P_nu;
}

void PMatrixBuilder::setupNormalizationConditions()
{
    m_normalization_conditions.clear();
    
    // Set up 3 normalization conditions per Newton step to fix degrees of freedom
    // Standard choices: fix one coefficient per boundary component
    
    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        if (nu == 0)
        {
            // First component: typically fix constant and linear terms
            m_normalization_conditions.push_back({nu, 0, 0.0});   // Constant term
            m_normalization_conditions.push_back({nu, 1, 1.0});   // Linear term
        }
        else if (nu == 1 && m_normalization_conditions.size() < 3)
        {
            // Second component: fix one more degree of freedom
            m_normalization_conditions.push_back({nu, -1, 0.0});  // First negative frequency
        }
    }
    
    // Debug logging removed for compilation
}

void PMatrixBuilder::validateParameters() const
{
    if (m_connectivity < 2)
    {
        throw std::invalid_argument("PMatrixBuilder: Connectivity must be at least 2");
    }
    
    if (m_N <= 0 || (m_N & (m_N - 1)) != 0)
    {
        throw std::invalid_argument("PMatrixBuilder: N must be a positive power of 2");
    }
    
    if (m_is_annulus && m_connectivity != 2)
    {
        throw std::invalid_argument("PMatrixBuilder: Annulus mode requires exactly 2 boundary components");
    }
}

int PMatrixBuilder::getPositiveFrequencyCount(int nu) const
{
    // Number of positive frequencies to zero for component nu
    // Typically this is N/2 - 1 (excluding zero frequency)
    return m_N/2 - 1;
}

int PMatrixBuilder::getNegativeFrequencyCount(int nu) const
{
    // Number of negative frequencies to keep for component nu
    // Typically this is N/2 for the Laurent series representation
    return m_N/2;
}