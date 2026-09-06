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

#pragma once

#include "../core/Types.h"
#include "MCSCReflection.h"

#include <vector>

class MCSCCircleDomain;
class MCSCPolygonalDomain;

/**
 * @brief MCSC reflection-method integrand for the unbounded (exterior) case (dissertation eq.
 *        2.13; port of fpextrefl.m).
 *
 * Builds the truncated product f'(z) that is later integrated (with vertex-singularity
 * quadrature, #164) to recover the map f(z). f'(z) is the actual "integrand" the algorithm
 * evaluates -- the singularity function S(z) = f''(z)/f'(z) discussed in the dissertation's
 * convergence proof (Sec 3.6) is never materialized separately, here or in the MATLAB reference.
 *
 * Unlike the MATLAB reference's addlistener-based auto-rebuild on circle-domain changes, this
 * class requires an explicit rebuild() call after the circle domain is mutated (e.g. by a Newton
 * step via MCSCCircleDomain::setFromUnconstrained) -- matching this project's existing
 * explicit-recompute convention (see FornbergMC).
 */
class MCSCReflectionIntegrand
{
public:
    /**
     * @brief Build the integrand for the given target/source domains and truncation level.
     * @param polygon Target (unbounded) polygonal domain -- supplies the vertex turning angles
     *        alpha used to compute beta(k,j) = 1 - alpha(k,j).
     * @param circle Source circle domain -- supplies circle centers/radii/prevertices.
     * @param N Reflection truncation level (N >= 0).
     * @throws std::invalid_argument if polygon and circle have different connectivity, or if any
     *         circle's prevertex count doesn't match its corresponding polygon's vertex count.
     */
    MCSCReflectionIntegrand(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N);

    /**
     * @brief Rebuild the reflection data after the circle domain has changed.
     *
     * Call this after any mutation of the circle domain used to construct this integrand (e.g.
     * MCSCCircleDomain::setFromUnconstrained during a Newton step) -- reflection data is not
     * automatically kept in sync.
     *
     * @throws std::invalid_argument under the same conditions as the constructor.
     */
    void rebuild(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N);

    /**
     * @brief Evaluate f'(z) (dissertation eq. 2.13; port of fpextrefl.eval_fprime).
     *
     * z coinciding exactly with a prevertex or one of its reflected images produces a
     * divide-by-zero (inf/nan result), matching the MATLAB reference's identical behavior --
     * not guarded here because the vertex-singularity quadrature (#164) is expected to
     * deliberately evaluate arbitrarily close to (though never exactly at) these points by
     * design, so a hard throw would work against its intended use rather than catch a bug.
     *
     * @param z Evaluation point (should not coincide with a prevertex or its reflected images).
     * @return f'(z).
     */
    Complex evalFPrime(const Complex& z) const;

    /**
     * @brief Evaluate the single-vertex partial factor f_{k,j}(z) (port of
     *        fpextrefl.eval_fkj) -- the contribution of vertex k on target component j alone,
     *        needed by the parameter-problem/continuation solver (#165).
     * @param z Evaluation point.
     * @param k Vertex index on component j (0-based).
     * @param j Boundary component index (0-based).
     * @return f_{k,j}(z).
     * @throws std::invalid_argument if j or k is out of range.
     */
    Complex evalFkj(const Complex& z, int k, int j) const;

private:
    std::vector<std::vector<double>> m_beta;                       // beta[j][k]
    std::vector<std::vector<mcsc::ReflectedCircle>> m_reflections;  // m_reflections[j][nu]

    static void validateDomains(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle);
};
