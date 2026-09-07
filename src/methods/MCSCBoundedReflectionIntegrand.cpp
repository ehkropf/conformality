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

#include "MCSCBoundedReflectionIntegrand.h"

#include "../domains/MCSCCircleDomain.h"
#include "../domains/MCSCPolygonalDomain.h"

#include <cmath>
#include <complex>
#include <stdexcept>

void MCSCBoundedReflectionIntegrand::validateDomains(
    const MCSCPolygonalDomain& polygon,
    const MCSCCircleDomain& circle
)
{
    if (polygon.getConnectivity() != circle.circleCount())
    {
        throw std::invalid_argument(
            "MCSCBoundedReflectionIntegrand: polygon and circle domain connectivity must match");
    }
    if (circle.circleCount() < 2)
    {
        throw std::invalid_argument("MCSCBoundedReflectionIntegrand: at least 2 circles are required");
    }

    for (int j = 0; j < circle.circleCount(); ++j)
    {
        if (polygon.vertexCount(j) != static_cast<int>(circle.prevertexCount(j)))
        {
            throw std::invalid_argument(
                "MCSCBoundedReflectionIntegrand: prevertex count on circle j must match vertex count on target component j");
        }
    }
}

MCSCBoundedReflectionIntegrand::MCSCBoundedReflectionIntegrand(
    const MCSCPolygonalDomain& polygon,
    const MCSCCircleDomain& circle,
    int N
)
{
    rebuild(polygon, circle, N);
}

void MCSCBoundedReflectionIntegrand::rebuild(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N)
{
    validateDomains(polygon, circle);

    const int m = circle.circleCount();

    // beta(k,j) = 1 - alpha(k,j), except the outer component (j=0, C1) whose sign is flipped
    // (intpolys.m: beta(:,1) = -beta(:,1), matching the dissertation's sum(beta_{k,1}) = -2
    // convention vs. +2 for every other component).
    std::vector<std::vector<double>> beta(m);
    for (int j = 0; j < m; ++j)
    {
        const auto& alpha = polygon.getAlpha(j);
        beta[j].resize(alpha.size());
        for (std::size_t k = 0; k < alpha.size(); ++k)
        {
            beta[j][k] = (j == 0) ? -(1.0 - alpha[k]) : (1.0 - alpha[k]);
        }
    }

    std::vector<Complex> centers = circle.getCenters();
    std::vector<double> radii = circle.getRadii();
    std::vector<std::vector<Complex>> prevertices(m);
    for (int j = 0; j < m; ++j)
    {
        prevertices[j] = circle.getPrevertices(j);
    }

    // Pre-reflection step (fpintrefl.build_reflections): reflect every interior circle (j >= 1)
    // through C1 (circle 0) to its exterior image before running the shared reflection
    // machinery. C1 itself is left untouched -- it is not reflected through itself.
    for (int j = 1; j < m; ++j)
    {
        const auto [co, ro] = mcsc::reflectCircle(centers[0], radii[0], centers[j], radii[j]);
        centers[j] = co;
        radii[j] = ro;
        for (Complex& z : prevertices[j])
        {
            z = mcsc::reflectPoint(centers[0], radii[0], z);
        }
    }

    // outerCenters mirrors the (post-pre-reflection) centers, matching the unbounded case's
    // outerImage = center convention -- there is no separate "s" quantity in the bounded case.
    std::vector<Complex> outerCenters = centers;

    m_beta = std::move(beta);
    m_reflections = mcsc::reflectCircleSequence(centers, radii, prevertices, outerCenters, N);
    m_N = N;
}

Complex MCSCBoundedReflectionIntegrand::evalFPrime(const Complex& z) const
{
    Complex logsum{0.0, 0.0};

    // Level-0 term: C1's own original, unreflected prevertices. m_reflections[0][0] is always
    // circle 0's untouched data (level 0 is never pre-reflected or otherwise altered).
    const auto& c1Level0 = m_reflections[0][0];
    for (std::size_t k = 0; k < c1Level0.prevertices.size(); ++k)
    {
        logsum += m_beta[0][k] * std::log(1.0 - (c1Level0.prevertices[k] - c1Level0.center) / (z - c1Level0.center));
    }

    // fpintrefl.m's eval_fprime "level" loop (1..N) is offset by one from reflectCircleSequence's
    // construction levels: its level=1 revisits construction-level-0 (each interior circle's own
    // pre-reflected-to-exterior copy, m_reflections[j][0] -- an ordinary-term contributor here,
    // unlike the initial draft of this class), level=2 covers construction-level-1's new entries,
    // and so on through construction-level (N-1). reflectCircleSequence's outermost,
    // construction-level-N entries are therefore built but never read -- a faithful port of
    // fpintrefl.m's own asymmetry between build_reflections's reflectzsmi(...,N) call and
    // eval_fprime's `for level=1:N`, not a bug to silently narrow away.
    //
    // nu1 is a single running counter walking C1's own reflection array (m_reflections[0]) by raw
    // sequential position -- NOT a matching reflection path -- one increment per (level, j, nu)
    // triple visited. Verified index-for-index against a direct simulation of reflectzsmi/
    // eval_fprime's index arithmetic for m=3,4,5 and N=1,2,3, and against the dissertation (Sec
    // 6.1 confirms the reflection method -- unlike the unresolved fast/Laurent-series method of
    // Ch. 4 -- was completed and validated for the bounded case).
    const int m = static_cast<int>(m_reflections.size());
    const auto& reflC1 = m_reflections[0];

    std::size_t levelStart = 0;
    std::size_t levelSize = 1;
    std::size_t nu1 = 0;
    for (int level = 0; level < m_N; ++level)
    {
        const std::size_t levelEnd = levelStart + levelSize;

        for (int j = 1; j < m; ++j)
        {
            const auto& reflJ = m_reflections[j];
            for (std::size_t nu = levelStart; nu < levelEnd; ++nu)
            {
                const auto& refl = reflJ[nu];
                for (std::size_t k = 0; k < refl.prevertices.size(); ++k)
                {
                    logsum += m_beta[j][k] * std::log(1.0 - (refl.prevertices[k] - refl.center) / (z - refl.center));
                }

                ++nu1;
                const auto& c1AtSamePosition = reflC1[nu1];
                for (std::size_t k = 0; k < c1AtSamePosition.prevertices.size(); ++k)
                {
                    logsum += m_beta[0][k]
                        * std::log(1.0 - (c1AtSamePosition.prevertices[k] - refl.center) / (z - refl.center));
                }
            }
        }

        levelStart = levelEnd;
        levelSize *= static_cast<std::size_t>(m - 1);
    }

    return std::exp(logsum);
}
