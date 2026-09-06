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

#include "MCSCUnconstrainedAngles.h"
#include "../core/Types.h"

#include <cmath>
#include <stdexcept>

namespace mcsc
{

std::vector<double> anglesToUnconstrained(const std::vector<double>& theta)
{
    const std::size_t K = theta.size();
    if (K < 2)
    {
        throw std::invalid_argument("anglesToUnconstrained: at least 2 prevertex angles are required");
    }

    // Gaps phi_k = theta_{k+1} - theta_k for k = 1, ..., K-1, and phi_K = (theta_1 + 2*pi) - theta_K,
    // matching circdomain.m's diff([theta(1); ...; theta(K)], theta(1)+2*pi).
    std::vector<double> phi(K);
    for (std::size_t k = 0; k + 1 < K; ++k)
    {
        phi[k] = theta[k + 1] - theta[k];
        if (phi[k] <= 0.0)
        {
            throw std::invalid_argument("anglesToUnconstrained: theta must be strictly increasing");
        }
    }
    phi[K - 1] = (theta[0] + TWO_PI) - theta[K - 1];
    if (phi[K - 1] <= 0.0)
    {
        throw std::invalid_argument("anglesToUnconstrained: theta must be strictly increasing and span less than 2*pi");
    }

    // psi_k = log(phi_{k+1} / phi_1) for k = 1, ..., K-1 (dissertation eq. 3.4).
    std::vector<double> psi(K - 1);
    for (std::size_t k = 0; k + 1 < K; ++k)
    {
        psi[k] = std::log(phi[k + 1] / phi[0]);
    }

    return psi;
}

std::vector<double> anglesFromUnconstrained(double theta1, const std::vector<double>& psi)
{
    if (psi.empty())
    {
        throw std::invalid_argument("anglesFromUnconstrained: at least 1 unconstrained variable is required");
    }

    const std::size_t K = psi.size() + 1;

    // sk(0) = 1, sk(k) = 1 + sum_{j=1}^{k} exp(psi_j), so sk(K-1) = 1 + sum(exp(psi)).
    // theta_k = theta1 + 2*pi * sk(k-2) / sk(K-2) for k = 2, ..., K  (dissertation eq. 3.5).
    std::vector<double> sk(K - 1);
    sk[0] = 1.0;
    for (std::size_t j = 0; j + 1 < K - 1; ++j)
    {
        sk[j + 1] = sk[j] + std::exp(psi[j]);
    }
    const double total = sk.back() + std::exp(psi.back());

    std::vector<double> theta(K);
    theta[0] = theta1;
    for (std::size_t k = 1; k < K; ++k)
    {
        theta[k] = theta1 + TWO_PI * sk[k - 1] / total;
    }

    return theta;
}

} // namespace mcsc
