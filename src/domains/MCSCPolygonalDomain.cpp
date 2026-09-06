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

#include "MCSCPolygonalDomain.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

MCSCPolygonalDomain::MCSCPolygonalDomain(
    const std::vector<std::vector<Complex>>& vertices,
    bool isUnboundedDomain
)
    : Domain(isUnboundedDomain, static_cast<int>(vertices.size()))
    , m_vertices{vertices}
{
    if (m_vertices.empty())
    {
        throw std::invalid_argument("MCSCPolygonalDomain: Must have at least one boundary component");
    }

    m_alpha.resize(m_vertices.size());
    for (size_t i = 0; i < m_vertices.size(); ++i)
    {
        if (m_vertices[i].size() < 3)
        {
            throw std::invalid_argument(
                "MCSCPolygonalDomain: Each boundary component must have at least 3 vertices");
        }

        m_alpha[i] = calcAngles(m_vertices[i]);
    }
}

const std::vector<Complex>& MCSCPolygonalDomain::getVertices(int componentIndex) const
{
    if (componentIndex < 0 || componentIndex >= static_cast<int>(m_vertices.size()))
    {
        throw std::invalid_argument("MCSCPolygonalDomain: Component index out of range");
    }
    return m_vertices[componentIndex];
}

const std::vector<double>& MCSCPolygonalDomain::getAlpha(int componentIndex) const
{
    if (componentIndex < 0 || componentIndex >= static_cast<int>(m_alpha.size()))
    {
        throw std::invalid_argument("MCSCPolygonalDomain: Component index out of range");
    }
    return m_alpha[componentIndex];
}

int MCSCPolygonalDomain::vertexCount(int componentIndex) const
{
    return static_cast<int>(getVertices(componentIndex).size());
}

bool MCSCPolygonalDomain::contains(const Complex& z) const
{
    const double boundaryTolerance = BOUNDARY_TOLERANCE;

    if (isUnbounded())
    {
        // Unbounded target domain: the region is the exterior of all m disjoint
        // polygons (no enclosing outer boundary), so z is inside the domain if and
        // only if it is outside every one of them.
        for (const auto& component : m_vertices)
        {
            int winding = Domain::calculateWindingNumber(z, component, boundaryTolerance);
            if (winding == std::numeric_limits<int>::max())
            {
                return false;
            }
            if (winding != 0)
            {
                return false;
            }
        }
        return true;
    }

    // Bounded target domain: component 0 is the enclosing outer boundary, and every
    // remaining component is a hole. z is inside the domain iff it is inside the
    // outer boundary and outside every hole.
    int outerWinding = Domain::calculateWindingNumber(z, m_vertices[0], boundaryTolerance);
    if (outerWinding == std::numeric_limits<int>::max())
    {
        return true;
    }
    if (outerWinding == 0)
    {
        return false;
    }

    for (size_t i = 1; i < m_vertices.size(); ++i)
    {
        int winding = Domain::calculateWindingNumber(z, m_vertices[i], boundaryTolerance);
        if (winding == std::numeric_limits<int>::max())
        {
            return false;
        }
        if (winding != 0)
        {
            return false;
        }
    }

    return true;
}

void MCSCPolygonalDomain::transformBoundary(std::function<Complex(const Complex&)> transform)
{
    // Compute into local copies and validate before committing, so a failure partway
    // through leaves this domain completely unmodified (strong exception safety).
    std::vector<std::vector<Complex>> newVertices = m_vertices;
    for (auto& component : newVertices)
    {
        for (auto& vertex : component)
        {
            vertex = transform(vertex);
        }
    }

    std::vector<std::vector<double>> newAlpha(newVertices.size());
    for (size_t i = 0; i < newVertices.size(); ++i)
    {
        newAlpha[i] = calcAngles(newVertices[i]);
    }

    m_vertices = std::move(newVertices);
    m_alpha = std::move(newAlpha);
}

std::vector<double> MCSCPolygonalDomain::calcAngles(std::vector<Complex>& vertices)
{
    const size_t n = vertices.size();

    std::vector<Complex> incoming(n);
    std::vector<Complex> outgoing(n);
    for (size_t j = 0; j < n; ++j)
    {
        const size_t prev = (j == 0) ? n - 1 : j - 1;
        incoming[j] = vertices[j] - vertices[prev];
    }
    for (size_t j = 0; j < n; ++j)
    {
        const size_t next = (j + 1) % n;
        outgoing[j] = incoming[next];
    }

    std::vector<double> alpha(n);
    for (size_t j = 0; j < n; ++j)
    {
        double turningAngle = std::arg(-incoming[j] * std::conj(outgoing[j]));
        double a = std::fmod(turningAngle / M_PI, 2.0);
        if (a < 0.0)
        {
            a += 2.0;
        }
        alpha[j] = a;
    }

    double angleSum = 0.0;
    for (double a : alpha)
    {
        angleSum += (1.0 - a);
    }

    const double tolerance = 100.0 * std::sqrt(static_cast<double>(n)) * std::numeric_limits<double>::epsilon();
    if (std::abs(angleSum - std::round(angleSum)) > tolerance)
    {
        throw std::invalid_argument("MCSCPolygonalDomain: Invalid polygon (interior angle sum check failed)");
    }

    if (angleSum < 0.0)
    {
        std::reverse(vertices.begin(), vertices.end());
        std::reverse(alpha.begin(), alpha.end());
    }

    return alpha;
}
