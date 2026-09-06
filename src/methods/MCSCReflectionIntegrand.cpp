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

#include "MCSCReflectionIntegrand.h"

#include "../domains/MCSCCircleDomain.h"
#include "../domains/MCSCPolygonalDomain.h"

#include <cmath>
#include <complex>
#include <stdexcept>

void MCSCReflectionIntegrand::validateDomains(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle)
{
    if (polygon.getConnectivity() != circle.circleCount())
    {
        throw std::invalid_argument(
            "MCSCReflectionIntegrand: polygon and circle domain connectivity must match");
    }

    for (int j = 0; j < circle.circleCount(); ++j)
    {
        if (polygon.vertexCount(j) != static_cast<int>(circle.prevertexCount(j)))
        {
            throw std::invalid_argument(
                "MCSCReflectionIntegrand: prevertex count on circle j must match vertex count on target component j");
        }
    }
}

MCSCReflectionIntegrand::MCSCReflectionIntegrand(
    const MCSCPolygonalDomain& polygon,
    const MCSCCircleDomain& circle,
    int N
)
{
    rebuild(polygon, circle, N);
}

void MCSCReflectionIntegrand::rebuild(const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& circle, int N)
{
    validateDomains(polygon, circle);

    const int m = circle.circleCount();

    // beta(k,j) = 1 - alpha(k,j), no outer-component sign flip for the unbounded case
    // (extpolys.m; contrast with intpolys.m's beta(:,1) = -beta(:,1) for the bounded case).
    std::vector<std::vector<double>> beta(m);
    for (int j = 0; j < m; ++j)
    {
        const auto& alpha = polygon.getAlpha(j);
        beta[j].resize(alpha.size());
        for (std::size_t k = 0; k < alpha.size(); ++k)
        {
            beta[j][k] = 1.0 - alpha[k];
        }
    }

    std::vector<Complex> centers = circle.getCenters();
    std::vector<double> radii = circle.getRadii();
    std::vector<std::vector<Complex>> prevertices(m);
    for (int j = 0; j < m; ++j)
    {
        prevertices[j] = circle.getPrevertices(j);
    }

    // Unbounded case (fpextrefl.build_reflections): the "outer image" reflected alongside the
    // prevertices is each circle's own center (cn passed as both cn and sn in the MATLAB code) --
    // no pre-reflection-to-exterior step, unlike the bounded case (#167).
    std::vector<Complex> outerCenters = centers;

    m_beta = std::move(beta);
    m_reflections = mcsc::reflectCircleSequence(centers, radii, prevertices, outerCenters, N);
}

Complex MCSCReflectionIntegrand::evalFPrime(const Complex& z) const
{
    Complex zprod1{1.0, 0.0};
    Complex logzprod2{0.0, 0.0};

    const int m = static_cast<int>(m_reflections.size());
    for (int j = 0; j < m; ++j)
    {
        for (const auto& refl : m_reflections[j])
        {
            const Complex zs = z - refl.center;
            for (std::size_t k = 0; k < refl.prevertices.size(); ++k)
            {
                logzprod2 += m_beta[j][k] * std::log(1.0 - (refl.prevertices[k] - refl.center) / zs);
            }
            // std::pow(complex, 2) goes through exp(2*log(z)), which is measurably less precise
            // than a direct square (e.g. introduces a spurious ~1e-15 imaginary part on a purely
            // real input) -- squaring directly avoids that avoidable error, accumulated over
            // every reflection node in the product.
            const Complex ratio = zs / (z - refl.outerImage);
            zprod1 *= ratio * ratio;
        }
    }

    return zprod1 * std::exp(logzprod2);
}

Complex MCSCReflectionIntegrand::evalFkj(const Complex& z, int k, int j) const
{
    const int m = static_cast<int>(m_reflections.size());
    if (j < 0 || j >= m)
    {
        throw std::invalid_argument("MCSCReflectionIntegrand::evalFkj: component index j out of range");
    }
    if (k < 0 || k >= static_cast<int>(m_beta[j].size()))
    {
        throw std::invalid_argument("MCSCReflectionIntegrand::evalFkj: vertex index k out of range");
    }

    Complex logprod{0.0, 0.0};
    for (const auto& refl : m_reflections[j])
    {
        logprod += std::log(1.0 - (refl.prevertices[static_cast<std::size_t>(k)] - refl.center) / (z - refl.center));
        logprod += std::log(1.0 - (refl.center - refl.outerImage) / (z - refl.outerImage));
    }

    return std::exp(logprod);
}
