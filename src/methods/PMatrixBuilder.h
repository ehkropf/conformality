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

#include "FornbergMCConfiguration.h"
#include <Eigen/Dense>
#include <vector>

/**
 * @brief Builder class for P_ν matrices in the Fornberg-like method
 *
 * This class constructs the P_ν matrices that enforce analyticity conditions
 * by zeroing positive frequency coefficients in the Laurent series representation.
 * It handles both general multiply connected cases and specialized annulus 
 * formulations for improved efficiency.
 *
 * The P_ν matrices are central to the Fornberg method's linear system formation,
 * ensuring that the computed boundary correspondences satisfy the required
 * analyticity conditions for conformal mappings.
 */
class PMatrixBuilder
{
private:
    const FornbergMCConfiguration& m_config;
    int m_connectivity;                    // Number of boundary components (m)
    int m_N;                              // Boundary points per component
    bool m_is_annulus;                    // Whether using annulus optimization
    
    // Cached matrices for efficiency
    std::vector<Eigen::MatrixXcd> m_P_matrices;  // P_ν matrices for each component
    std::vector<Eigen::VectorXi> m_frequency_indices; // Frequency indexing per component
    
    // Normalization information
    struct NormalizationCondition
    {
        int component;      // Which boundary component
        int frequency;      // Which frequency to normalize
        double value;       // Normalization value
    };
    std::vector<NormalizationCondition> m_normalization_conditions;

public:
    /**
     * @brief Construct a new P Matrix Builder
     * @param config Configuration parameters
     * @param connectivity Number of boundary components
     * @param is_annulus Whether to use annulus optimization
     */
    PMatrixBuilder(const FornbergMCConfiguration& config, int connectivity, bool is_annulus = false);

    /**
     * @brief Destructor
     */
    ~PMatrixBuilder() = default;

    /**
     * @brief Build P_ν matrix for a specific boundary component
     * @param nu Component index (0-based)
     * @return P_ν matrix for component nu
     */
    Eigen::MatrixXcd buildPMatrix(int nu) const;

    /**
     * @brief Build all P_ν matrices
     * @return Vector of P_ν matrices for all components
     */
    std::vector<Eigen::MatrixXcd> buildAllPMatrices() const;

    /**
     * @brief Get frequency indices for a boundary component
     * @param nu Component index (0-based)
     * @return Vector of frequency indices used for component nu
     */
    const Eigen::VectorXi& getFrequencyIndices(int nu) const;

    /**
     * @brief Apply normalization conditions to the system
     * @param system_matrix System matrix to modify
     * @param rhs_vector Right-hand side vector to modify
     */
    void applyNormalizationConditions(Eigen::MatrixXcd& system_matrix, 
                                    Eigen::VectorXcd& rhs_vector) const;

    /**
     * @brief Get the size of the constructed system
     * @return Total number of equations in the P-matrix system
     */
    int getSystemSize() const;

    /**
     * @brief Check if using annulus optimization
     * @return True if annulus formulation is active
     */
    bool isAnnulusMode() const
    {
        return m_is_annulus;
    }

    /**
     * @brief Update connectivity (requires rebuilding matrices)
     * @param new_connectivity New number of boundary components
     */
    void setConnectivity(int new_connectivity);

    /**
     * @brief Switch between general and annulus modes
     * @param use_annulus Whether to use annulus optimization
     */
    void setAnnulusMode(bool use_annulus);

private:
    /**
     * @brief Initialize frequency indexing for all components
     */
    void initializeFrequencyIndices();

    /**
     * @brief Build P matrix for general multiply connected case
     * @param nu Component index
     * @return P_ν matrix for general case
     */
    Eigen::MatrixXcd buildGeneralPMatrix(int nu) const;

    /**
     * @brief Build P matrix for annulus case (optimized)
     * @param nu Component index (should be 0 or 1 for annulus)
     * @return P_ν matrix for annulus case
     */
    Eigen::MatrixXcd buildAnnulusPMatrix(int nu) const;

    /**
     * @brief Set up normalization conditions based on connectivity
     */
    void setupNormalizationConditions();

    /**
     * @brief Validate matrix construction parameters
     * @throws std::invalid_argument if parameters are inconsistent
     */
    void validateParameters() const;

    /**
     * @brief Get the number of positive frequencies to zero for component nu
     * @param nu Component index
     * @return Number of positive frequencies to enforce as zero
     */
    int getPositiveFrequencyCount(int nu) const;

    /**
     * @brief Get the number of negative frequencies to keep for component nu
     * @param nu Component index  
     * @return Number of negative frequencies to retain
     */
    int getNegativeFrequencyCount(int nu) const;
};