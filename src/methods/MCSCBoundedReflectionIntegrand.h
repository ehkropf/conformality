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
 * @brief MCSC reflection-method integrand for the bounded (interior) case (dissertation eq.
 *        2.16/3.2; port of fpintrefl.m).
 *
 * Builds the truncated product f'(z) that is later integrated (with vertex-singularity
 * quadrature, #164) to recover the map f(z). Unlike the unbounded case's
 * MCSCReflectionIntegrand, there is no s_nu-denominator anywhere -- the bounded product is a
 * bare polynomial-type factor -- and circle 0 (C1, the outer circle) requires an "extra level of
 * reflection" pairing (dissertation remark after eq. 2.18/3.2) to make the truncated singularity
 * function's residues cancel in pairs.
 *
 * Like MCSCReflectionIntegrand, this class requires an explicit rebuild() call after the circle
 * domain is mutated (e.g. by a Newton step) -- matching this project's existing
 * explicit-recompute convention (see FornbergMC) rather than the MATLAB reference's
 * addlistener-based auto-rebuild.
 */
class MCSCBoundedReflectionIntegrand
{
public:
    /**
     * @brief Build the integrand for the given target/source domains and truncation level.
     * @param polygon Target (bounded) polygonal domain -- supplies the vertex turning angles
     *        alpha used to compute beta(k,j). beta(k,0) = -(1 - alpha(k,0)) for the outer
     *        component (sign-flipped, per intpolys.m); beta(k,j) = 1 - alpha(k,j) for j >= 1.
     * @param circle Source circle domain -- supplies circle centers/radii/prevertices. Circle 0
     *        is the outer circle C1; circles 1, ..., m-1 are interior circles.
     * @param N Reflection truncation level (N >= 0).
     * @throws std::invalid_argument if polygon and circle have different connectivity, if any
     *         circle's prevertex count doesn't match its corresponding polygon's vertex count,
     *         if there are fewer than 2 circles, or if N is negative (propagated from
     *         mcsc::reflectCircleSequence).
     */
    MCSCBoundedReflectionIntegrand(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N);

    /**
     * @brief Rebuild the reflection data after the circle domain has changed.
     *
     * Call this after any mutation of the circle domain used to construct this integrand (e.g.
     * MCSCCircleDomain::setFromUnconstrained during a Newton step) -- reflection data is not
     * automatically kept in sync.
     *
     * @throws std::invalid_argument under the same conditions as the constructor. On throw, this
     *         object is left unchanged (strong exception guarantee) -- m_beta, m_reflections, and
     *         m_N are only overwritten together, after the throwing steps have succeeded.
     */
    void rebuild(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N);

    /**
     * @brief Evaluate f'(z) (dissertation eq. 2.16/3.2; port of fpintrefl.eval_fprime).
     *
     * z coinciding exactly with a prevertex or one of its reflected images produces a
     * divide-by-zero (inf/nan result), matching the MATLAB reference's identical behavior --
     * not guarded here for the same reason as MCSCReflectionIntegrand::evalFPrime: the
     * vertex-singularity quadrature (#164) is expected to deliberately evaluate arbitrarily
     * close to (though never exactly at) these points by design.
     *
     * @param z Evaluation point (should not coincide with a prevertex or its reflected images).
     * @return f'(z).
     */
    Complex evalFPrime(const Complex& z) const;

private:
    std::vector<std::vector<double>> m_beta;                       // beta[j][k]
    std::vector<std::vector<mcsc::ReflectedCircle>> m_reflections;  // m_reflections[j][nu]

    // Truncation level bounding evalFPrime's level loop -- NOT redundant with m_reflections'
    // size: m_reflections is built one construction level deeper than evalFPrime ever reads (see
    // evalFPrime's comment), so this is load-bearing, separate state. Set together with
    // m_reflections in rebuild() only; must never be modified independently.
    int m_N{0};

    static void validateDomains(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle);
};
