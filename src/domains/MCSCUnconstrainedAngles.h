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

#include <vector>

/**
 * @brief Unconstrained-parameter transform for ordered prevertex angles (dissertation Sec 3.2).
 *
 * The prevertex angles theta_{1,j} < theta_{2,j} < ... < theta_{Kj,j} on a single circle must
 * satisfy the wraparound sum constraint sum_k (theta_{k+1} - theta_k) = 2*pi (with the convention
 * theta_{K+1} = theta_1 + 2*pi). This ordering/sum constraint is difficult to enforce directly in
 * a Newton solve, so MCSC (following [Trefethen, Driscoll et al.]) transforms the gaps
 * phi_k := theta_{k+1} - theta_k into unconstrained log-ratio variables
 *
 *   psi_k := log(phi_{k+1} / phi_1),  k = 1, ..., K-1
 *
 * which can take any real value while phi_k = phi_1 * exp(psi_{k-1}) > 0 always holds, so the
 * angles recovered from psi are automatically strictly increasing and sum to 2*pi. This mirrors
 * circdomain.m's unconstrained()/set.Xu() angle handling.
 *
 * These free functions operate on a single circle's angles in isolation (theta_1 assumed known /
 * fixed by the caller, exactly as in the dissertation and circdomain.m -- on circle 1, theta_1 is
 * fixed at 0; on other circles theta_1,j itself becomes an additional unconstrained parameter that
 * MCSCCircleDomain packs separately). They are intentionally decoupled from MCSCCircleDomain so
 * the transform itself is independently unit-testable.
 */
namespace mcsc
{

/**
 * @brief Transform K ordered prevertex angles theta_1 < theta_2 < ... < theta_K < theta_1 + 2*pi
 *        into K-1 unconstrained log-ratio variables psi_1, ..., psi_{K-1}.
 * @param theta Ordered prevertex angles for one circle (size K >= 2), theta_1 fixed by convention.
 * @return Unconstrained variables psi (size K-1).
 * @throws std::invalid_argument if theta has fewer than 2 elements or is not strictly increasing.
 */
std::vector<double> anglesToUnconstrained(const std::vector<double>& theta);

/**
 * @brief Inverse of anglesToUnconstrained(): recover ordered angles from theta_1 and psi.
 * @param theta1 First prevertex angle theta_1 (fixed / known).
 * @param psi Unconstrained log-ratio variables (size K-1).
 * @return Ordered prevertex angles theta_1, ..., theta_K (size K).
 */
std::vector<double> anglesFromUnconstrained(double theta1, const std::vector<double>& psi);

} // namespace mcsc
