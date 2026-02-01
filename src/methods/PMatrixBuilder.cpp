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

Eigen::MatrixXcd PMatrixBuilder::buildPMatrix(int nu, const ConformalModuli& moduli) const
{
    if (nu < 0 || nu >= m_connectivity)
    {
        throw std::invalid_argument("PMatrixBuilder: Invalid component index " + std::to_string(nu));
    }

    if (m_is_annulus)
    {
        return buildAnnulusPMatrix(nu, moduli);
    }
    else
    {
        return buildGeneralPMatrix(nu, moduli);
    }
}

std::vector<Eigen::MatrixXcd> PMatrixBuilder::buildAllPMatrices(const ConformalModuli& moduli) const
{
    std::vector<Eigen::MatrixXcd> P_matrices;
    P_matrices.reserve(m_connectivity);

    for (int nu = 0; nu < m_connectivity; ++nu)
    {
        P_matrices.push_back(buildPMatrix(nu, moduli));
    }

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
                                                  Eigen::VectorXcd& rhs_vector,
                                                  double norm_cond_value) const
{
    // Annulus case: no extra normalization rows needed
    // The c(1)=0 constraint already fixes degrees of freedom
    if (m_is_annulus)
    {
        return;
    }

    // General case: add final normalization row
    // D(end,:) = [1, 0, 0, ..., 0]
    int last_row = system_matrix.rows() - 1;
    system_matrix.row(last_row).setZero();
    system_matrix(last_row, 0) = 1.0;

    // Adjust RHS (matches MATLAB namap.m lines 109-111):
    // g(end-1) += N * norm_cond(3)
    // g(end) = 0 (from [g_; 0])
    rhs_vector(last_row - 1) += static_cast<double>(m_N) * norm_cond_value;
    rhs_vector(last_row) = 0.0;
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
            if (j < m_N/2)
            {
                m_frequency_indices[nu][j] = j;  // Positive frequencies and zero
            }
            else
            {
                m_frequency_indices[nu][j] = j - m_N;  // Negative frequencies (including Nyquist)
            }
        }
    }
    
}

void PMatrixBuilder::doTri(Eigen::MatrixXcd& P_, double rho) const
{
    // Implements MATLAB make_Pnu.m do_tri helper (lines 56-63)
    // Expects first column of P_ to be pre-filled
    int M = P_.rows();

    // Column 1 (0-based): P_(:,2) = [0; rho*(1:M-1)'.*P_(1:M-1,1)]
    // Row 0 is explicitly 0, rows 1 to M-1 get filled
    P_(0, 1) = 0.0;
    for (int i = 1; i < M; ++i)
    {
        // MATLAB: rho * k * P_(k,1) where k = 1:M-1
        // C++ (0-based): rho * i * P_(i-1, 0) where i = 1:M-1
        P_(i, 1) = rho * static_cast<double>(i) * P_(i - 1, 0);
    }

    // Columns 2 to M-1 (0-based): j = 3:M in MATLAB becomes j = 2:M-1 in C++
    for (int j = 2; j < M; ++j)
    {
        // P_(j:M,j) = P_(j-1:M-1,j-1).*(j-1:M-1)'*(rho/(j-1))
        // In C++ (0-based): rows j to M-1, column j
        for (int i = j; i < M; ++i)
        {
            // MATLAB: P_(k,j) = P_(k-1,j-1) * (k-1) * rho / (j-1)
            // In MATLAB k goes j:M, so k-1 goes (j-1):(M-1)
            // In C++ (0-based), i goes j:M-1, and the multiplier is i (since MATLAB k = i+1, k-1 = i)
            P_(i, j) = P_(i - 1, j - 1) * static_cast<double>(i) * (rho / static_cast<double>(j));
        }
    }
}

Eigen::MatrixXcd PMatrixBuilder::buildGeneralPMatrix(int nu, const ConformalModuli& moduli) const
{
    // Implements MATLAB make_Pnu.m
    int M = m_N / 2;
    Eigen::MatrixXcd P_nu = Eigen::MatrixXcd::Zero(m_connectivity * M, m_N);

    if (nu == 0)
    {
        // First boundary component (outer boundary)
        // MATLAB: nu==1 case

        // First block-row: identity in upper-right block
        // MATLAB: Pnu(1:M, M+1:N) = eye(M)
        P_nu.block(0, M, M, M) = Eigen::MatrixXcd::Identity(M, M);

        // Rest of the block-rows
        // MATLAB: for L=1:conn-1
        for (int L = 0; L < m_connectivity - 1; ++L)
        {
            Eigen::MatrixXcd P_ = Eigen::MatrixXcd::Zero(M, M);

            // MATLAB: P_(:,1) = c(L).^(0:M-1)
            std::complex<double> c_L = moduli.c(L);
            for (int i = 0; i < M; ++i)
            {
                P_(i, 0) = std::pow(c_L, i);
            }

            // MATLAB: P_ = do_tri(P_, rho(L))
            doTri(P_, moduli.rho(L));

            // MATLAB: Pnu(L*M+1:(L+1)*M, 1:M) = P_.'
            // C++ (0-based): block((L+1)*M, 0, M, M) = P_.transpose()
            P_nu.block((L + 1) * M, 0, M, M) = P_.transpose();
        }
    }
    else
    {
        // Inner boundary components (nu >= 1 in C++, nu >= 2 in MATLAB)

        // First block-row: -rot90(P_, 2) in upper-right block
        Eigen::MatrixXcd P_ = Eigen::MatrixXcd::Zero(M, M);

        // MATLAB: P_(:,1) = rho(nu-1) * c(nu-1).^(0:M-1)
        std::complex<double> c_nu = moduli.c(nu - 1);
        double rho_nu = moduli.rho(nu - 1);
        for (int i = 0; i < M; ++i)
        {
            P_(i, 0) = rho_nu * std::pow(c_nu, i);
        }

        // MATLAB: P_ = do_tri(P_, rho(nu-1))
        doTri(P_, rho_nu);

        // MATLAB: Pnu(1:M, M+1:N) = -rot90(P_, 2)
        // rot90(P_, 2) reverses both rows and columns (180 degree rotation)
        P_nu.block(0, M, M, M) = -(P_.colwise().reverse().rowwise().reverse());

        // MATLAB: Pnu((nu-1)*M+1:nu*M, 1:M) = -eye(M)
        // C++ (0-based): block(nu*M, 0, M, M) = -Identity
        P_nu.block(nu * M, 0, M, M) = -Eigen::MatrixXcd::Identity(M, M);

        // Interaction blocks for other inner boundaries
        // MATLAB: for L=1:conn-1
        for (int L = 0; L < m_connectivity - 1; ++L)
        {
            // MATLAB: if L==(nu-1), continue, end
            if (L == nu - 1)
            {
                continue;
            }

            // MATLAB interaction formula:
            // pl_cvl = rho(L) / (c(nu-1) - c(L))
            // P_(1,:) = (1/(c(nu-1)-c(L))) * (rho(nu-1).^(1:M) ./ (c(L)-c(nu-1)).^(0:M-1))
            // P_(2,:) = pl_cvl * (1:M) .* P_(1,:)
            // for k=3:M, for j=1:M: P_(k,j) = pl_cvl * P_(k-1,j) * (j+k-2)/(k-1)

            std::complex<double> c_L = moduli.c(L);
            double rho_L = moduli.rho(L);
            std::complex<double> c_diff = c_nu - c_L;

            // Validate that hole centers are not coincident
            if (std::abs(c_diff) < 1e-14)
            {
                throw std::invalid_argument(
                    "PMatrixBuilder::buildGeneralPMatrix: Inner boundary centers c(" + std::to_string(nu - 1) +
                    ") and c(" + std::to_string(L) + ") are too close (distance " +
                    std::to_string(std::abs(c_diff)) + " < 1e-14). " +
                    "Ensure circle centers are distinct. Domain may be degenerate.");
            }
            std::complex<double> pl_cvl = rho_L / c_diff;

            Eigen::MatrixXcd P_int = Eigen::MatrixXcd::Zero(M, M);

            // First row: P_(1,:) = (1/(c_nu - c_L)) * (rho_nu^(1:M) / (c_L - c_nu)^(0:M-1))
            std::complex<double> neg_c_diff = c_L - c_nu;
            for (int j = 0; j < M; ++j)
            {
                // MATLAB: rho(nu-1)^(j+1) / (c(L)-c(nu-1))^j  where j is 0-based
                P_int(0, j) = (1.0 / c_diff) * std::pow(rho_nu, j + 1) / std::pow(neg_c_diff, j);
            }

            // Second row: P_(2,:) = pl_cvl * (1:M) .* P_(1,:)
            for (int j = 0; j < M; ++j)
            {
                // MATLAB: pl_cvl * (j+1) * P_(1, j+1) where j is 0-based
                P_int(1, j) = pl_cvl * static_cast<double>(j + 1) * P_int(0, j);
            }

            // Remaining rows: recurrence relation
            // MATLAB: for k=3:M, for j=1:M: P_(k,j) = pl_cvl * P_(k-1,j) * (j+k-2)/(k-1)
            for (int k = 2; k < M; ++k)  // k=2:M-1 in C++ (0-based) = k=3:M in MATLAB (1-based)
            {
                for (int j = 0; j < M; ++j)
                {
                    // MATLAB k and j are 1-based, so (j+k-2) with 1-based indices
                    // In C++ with 0-based: (j+1) + (k+1) - 2 = j + k
                    // MATLAB (k-1) with 1-based becomes k in C++ (0-based)
                    P_int(k, j) = pl_cvl * P_int(k - 1, j) * static_cast<double>(j + k) / static_cast<double>(k);
                }
            }

            // MATLAB: Pnu(L*M+1:(L+1)*M, M+1:N) = -fliplr(P_)
            // C++ (0-based): block((L+1)*M, M, M, M) = -P_.rowwise().reverse()
            P_nu.block((L + 1) * M, M, M, M) = -(P_int.rowwise().reverse());
        }
    }

    return P_nu;
}

Eigen::MatrixXcd PMatrixBuilder::buildAnnulusPMatrix(int nu, const ConformalModuli& moduli) const
{
    // Implements MATLAB make_Pnu_ann.m
    // Optimized structure for annulus case (m=2) and extends to higher connectivity
    int M = m_N / 2;
    Eigen::MatrixXcd P_nu = Eigen::MatrixXcd::Zero(m_connectivity * M, m_N);

    if (nu == 0)
    {
        // First boundary component (outer boundary)
        // MATLAB: nu==1 case

        // First block-row: identity (same as general)
        // MATLAB: Pnu(1:M, M+1:N) = eye(M)
        P_nu.block(0, M, M, M) = Eigen::MatrixXcd::Identity(M, M);

        // Second block-row: ANNULUS OPTIMIZATION
        // MATLAB: Pnu(M+1:2*M, 1:M) = diag(rho(1).^(0:M-1))
        double rho_1 = moduli.rho(0);
        for (int k = 0; k < M; ++k)
        {
            P_nu(M + k, k) = std::pow(rho_1, k);
        }

        // Rest of block-rows (for conn > 2, same as general)
        // MATLAB: for L=2:conn-1
        for (int L = 1; L < m_connectivity - 1; ++L)
        {
            Eigen::MatrixXcd P_ = Eigen::MatrixXcd::Zero(M, M);

            // MATLAB: P_(:,1) = c(L).^(0:M-1)
            std::complex<double> c_L = moduli.c(L);
            for (int i = 0; i < M; ++i)
            {
                P_(i, 0) = std::pow(c_L, i);
            }

            // MATLAB: P_ = do_tri(P_, rho(L))
            doTri(P_, moduli.rho(L));

            // MATLAB: Pnu(L*M+1:(L+1)*M, 1:M) = P_.'
            P_nu.block((L + 1) * M, 0, M, M) = P_.transpose();
        }
    }
    else if (nu == 1)
    {
        // First inner boundary (ANNULUS OPTIMIZATION for nu==2 in MATLAB)

        // First block-row: ANNULUS OPTIMIZATION
        // MATLAB: P_ = diag(rho(1).^(M:-1:1))
        // MATLAB: Pnu(1:M, M+1:N) = -P_
        double rho_1 = moduli.rho(0);
        for (int k = 0; k < M; ++k)
        {
            // MATLAB: rho(1)^(M:-1:1) puts rho^M at position 1, rho^1 at position M
            // C++ (0-based): row k gets -rho^(M-k)
            P_nu(k, M + k) = -std::pow(rho_1, M - k);
        }

        // Identity block (same as general)
        // MATLAB: Pnu((nu-1)*M+1:nu*M, 1:M) = -eye(M)
        P_nu.block(M, 0, M, M) = -Eigen::MatrixXcd::Identity(M, M);

        // Rest of block-rows (for conn > 2)
        // MATLAB: for L=2:conn-1
        for (int L = 1; L < m_connectivity - 1; ++L)
        {
            // Skip L=nu-1 (same as general)
            if (L == nu - 1)
            {
                continue;
            }

            // Interaction formula (same as general case)
            std::complex<double> c_nu = moduli.c(nu - 1);
            double rho_nu = moduli.rho(nu - 1);
            std::complex<double> c_L = moduli.c(L);
            double rho_L = moduli.rho(L);
            std::complex<double> c_diff = c_nu - c_L;

            // Validate that hole centers are not coincident
            if (std::abs(c_diff) < 1e-14)
            {
                throw std::invalid_argument(
                    "PMatrixBuilder::buildAnnulusPMatrix: Inner boundary centers c(" + std::to_string(nu - 1) +
                    ") and c(" + std::to_string(L) + ") are too close (distance " +
                    std::to_string(std::abs(c_diff)) + " < 1e-14). " +
                    "Ensure circle centers are distinct. Domain may be degenerate.");
            }
            std::complex<double> pl_cvl = rho_L / c_diff;

            Eigen::MatrixXcd P_int = Eigen::MatrixXcd::Zero(M, M);

            std::complex<double> neg_c_diff = c_L - c_nu;
            for (int j = 0; j < M; ++j)
            {
                P_int(0, j) = (1.0 / c_diff) * std::pow(rho_nu, j + 1) / std::pow(neg_c_diff, j);
            }

            for (int j = 0; j < M; ++j)
            {
                P_int(1, j) = pl_cvl * static_cast<double>(j + 1) * P_int(0, j);
            }

            for (int k = 2; k < M; ++k)
            {
                for (int j = 0; j < M; ++j)
                {
                    P_int(k, j) = pl_cvl * P_int(k - 1, j) * static_cast<double>(j + k) / static_cast<double>(k);
                }
            }

            P_nu.block((L + 1) * M, M, M, M) = -(P_int.rowwise().reverse());
        }
    }
    else
    {
        // Higher inner boundaries (nu >= 2 in C++, nu >= 3 in MATLAB)
        // These use the general structure but with c(2)=0 assumed real

        // First block-row (standard, not optimized)
        Eigen::MatrixXcd P_ = Eigen::MatrixXcd::Zero(M, M);
        std::complex<double> c_nu = moduli.c(nu - 1);
        double rho_nu = moduli.rho(nu - 1);

        for (int i = 0; i < M; ++i)
        {
            P_(i, 0) = rho_nu * std::pow(c_nu, i);
        }
        doTri(P_, rho_nu);

        // MATLAB: -fliplr(flipud(P_)) which is equivalent to -rot90(P_, 2)
        P_nu.block(0, M, M, M) = -(P_.colwise().reverse().rowwise().reverse());

        // Identity block
        P_nu.block(nu * M, 0, M, M) = -Eigen::MatrixXcd::Identity(M, M);

        // Interaction blocks
        // MATLAB: for L=2:conn-1 (note: starts at L=2 in annulus case)
        for (int L = 1; L < m_connectivity - 1; ++L)
        {
            if (L == nu - 1)
            {
                continue;
            }

            std::complex<double> c_L = moduli.c(L);
            double rho_L = moduli.rho(L);
            std::complex<double> c_diff = c_nu - c_L;

            // Validate that hole centers are not coincident (higher inner boundaries)
            if (std::abs(c_diff) < 1e-14)
            {
                throw std::invalid_argument(
                    "PMatrixBuilder::buildAnnulusPMatrix: Inner boundary centers c(" + std::to_string(nu - 1) +
                    ") and c(" + std::to_string(L) + ") are too close (distance " +
                    std::to_string(std::abs(c_diff)) + " < 1e-14). " +
                    "Ensure circle centers are distinct. Domain may be degenerate.");
            }
            std::complex<double> pl_cvl = rho_L / c_diff;

            Eigen::MatrixXcd P_int = Eigen::MatrixXcd::Zero(M, M);

            std::complex<double> neg_c_diff = c_L - c_nu;
            for (int j = 0; j < M; ++j)
            {
                P_int(0, j) = (1.0 / c_diff) * std::pow(rho_nu, j + 1) / std::pow(neg_c_diff, j);
            }

            for (int j = 0; j < M; ++j)
            {
                P_int(1, j) = pl_cvl * static_cast<double>(j + 1) * P_int(0, j);
            }

            for (int k = 2; k < M; ++k)
            {
                for (int j = 0; j < M; ++j)
                {
                    P_int(k, j) = pl_cvl * P_int(k - 1, j) * static_cast<double>(j + k) / static_cast<double>(k);
                }
            }

            P_nu.block((L + 1) * M, M, M, M) = -(P_int.rowwise().reverse());
        }
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
