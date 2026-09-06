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

#include "MCSCCircleDomain.h"
#include "MCSCUnconstrainedAngles.h"

#include <cmath>
#include <stdexcept>

namespace
{
constexpr double TWO_PI = 2.0 * M_PI;
}

MCSCCircleDomain::MCSCCircleDomain(std::vector<CircleData> circles)
    : Domain(false, static_cast<int>(circles.size()))
    , m_circles{std::move(circles)}
{
    if (m_circles.empty())
    {
        throw std::invalid_argument("MCSCCircleDomain: at least one circle is required");
    }

    for (const auto& circle : m_circles)
    {
        validateCircleData(circle);
    }
}

void MCSCCircleDomain::validateCircleData(const CircleData& circle)
{
    if (circle.radius <= 0.0)
    {
        throw std::invalid_argument("MCSCCircleDomain: circle radius must be positive");
    }

    const auto& t = circle.prevertexAngles;
    if (t.size() < 2)
    {
        throw std::invalid_argument("MCSCCircleDomain: circle must have at least 2 prevertices");
    }

    for (std::size_t k = 0; k + 1 < t.size(); ++k)
    {
        if (t[k + 1] <= t[k])
        {
            throw std::invalid_argument("MCSCCircleDomain: prevertex angles must be strictly increasing");
        }
    }
    if (t.back() >= t.front() + TWO_PI)
    {
        throw std::invalid_argument("MCSCCircleDomain: prevertex angles must span less than 2*pi");
    }
}

void MCSCCircleDomain::validateIndex(int j) const
{
    if (j < 0 || j >= circleCount())
    {
        throw std::invalid_argument("MCSCCircleDomain: circle index out of range");
    }
}

std::size_t MCSCCircleDomain::prevertexCount(int j) const
{
    validateIndex(j);
    return m_circles[j].prevertexAngles.size();
}

Complex MCSCCircleDomain::getCenter(int j) const
{
    validateIndex(j);
    return m_circles[j].center;
}

std::vector<Complex> MCSCCircleDomain::getCenters() const
{
    std::vector<Complex> centers;
    centers.reserve(m_circles.size());
    for (const auto& circle : m_circles)
    {
        centers.push_back(circle.center);
    }
    return centers;
}

double MCSCCircleDomain::getRadius(int j) const
{
    validateIndex(j);
    return m_circles[j].radius;
}

std::vector<double> MCSCCircleDomain::getRadii() const
{
    std::vector<double> radii;
    radii.reserve(m_circles.size());
    for (const auto& circle : m_circles)
    {
        radii.push_back(circle.radius);
    }
    return radii;
}

const std::vector<double>& MCSCCircleDomain::getPrevertexAngles(int j) const
{
    validateIndex(j);
    return m_circles[j].prevertexAngles;
}

std::vector<Complex> MCSCCircleDomain::getPrevertices(int j) const
{
    validateIndex(j);
    const auto& circle = m_circles[j];

    std::vector<Complex> prevertices;
    prevertices.reserve(circle.prevertexAngles.size());
    for (double t : circle.prevertexAngles)
    {
        prevertices.push_back(circle.center + circle.radius * Complex(std::cos(t), std::sin(t)));
    }
    return prevertices;
}

void MCSCCircleDomain::setCenter(int j, const Complex& center)
{
    validateIndex(j);
    m_circles[j].center = center;
}

void MCSCCircleDomain::setRadius(int j, double radius)
{
    validateIndex(j);
    if (radius <= 0.0)
    {
        throw std::invalid_argument("MCSCCircleDomain: circle radius must be positive");
    }
    m_circles[j].radius = radius;
}

void MCSCCircleDomain::setPrevertexAngles(int j, std::vector<double> angles)
{
    validateIndex(j);
    if (angles.size() != m_circles[j].prevertexAngles.size())
    {
        throw std::invalid_argument("MCSCCircleDomain: cannot change the number of prevertices on a circle");
    }

    CircleData candidate = m_circles[j];
    candidate.prevertexAngles = std::move(angles);
    validateCircleData(candidate);

    m_circles[j] = std::move(candidate);
}

int MCSCCircleDomain::expectedUnconstrainedSize() const
{
    const int m = circleCount();
    int size = 3 * m - 4;
    for (const auto& circle : m_circles)
    {
        size += static_cast<int>(circle.prevertexAngles.size());
    }
    return size;
}

Eigen::VectorXd MCSCCircleDomain::toUnconstrained() const
{
    const int m = circleCount();
    Eigen::VectorXd Xu(expectedUnconstrainedSize());

    // 1. log-radii of circles 1..m-1.
    for (int j = 1; j < m; ++j)
    {
        Xu[j - 1] = std::log(m_circles[j].radius);
    }

    // 2. Interleaved Re/Im of centers of circles 1..m-1.
    for (int j = 1; j < m; ++j)
    {
        const int base = (m - 1) + 2 * (j - 1);
        Xu[base] = std::real(m_circles[j].center);
        Xu[base + 1] = std::imag(m_circles[j].center);
    }

    // 3. Unconstrained angle variables for circle 0 (first prevertex fixed at its own angle).
    int b = 3 * m - 3;
    {
        const auto psi = mcsc::anglesToUnconstrained(m_circles[0].prevertexAngles);
        for (std::size_t k = 0; k < psi.size(); ++k)
        {
            Xu[b + static_cast<int>(k)] = psi[k];
        }
        b += static_cast<int>(psi.size());
    }

    // 4. For circles 1..m-1: first prevertex angle, then its unconstrained angle variables.
    for (int j = 1; j < m; ++j)
    {
        const auto& t = m_circles[j].prevertexAngles;
        Xu[b] = t.front();
        const auto psi = mcsc::anglesToUnconstrained(t);
        for (std::size_t k = 0; k < psi.size(); ++k)
        {
            Xu[b + 1 + static_cast<int>(k)] = psi[k];
        }
        b += 1 + static_cast<int>(psi.size());
    }

    return Xu;
}

void MCSCCircleDomain::setFromUnconstrained(const Eigen::VectorXd& Xu)
{
    const int m = circleCount();
    if (static_cast<int>(Xu.size()) != expectedUnconstrainedSize())
    {
        throw std::invalid_argument("MCSCCircleDomain: unconstrained parameter vector has wrong length");
    }

    // 1. Radii of circles 1..m-1 (circle 0's radius stays fixed at its current value, 1 by
    //    convention).
    for (int j = 1; j < m; ++j)
    {
        m_circles[j].radius = std::exp(Xu[j - 1]);
    }

    // 2. Centers of circles 1..m-1 (circle 0's center stays fixed by convention).
    for (int j = 1; j < m; ++j)
    {
        const int base = (m - 1) + 2 * (j - 1);
        m_circles[j].center = Complex(Xu[base], Xu[base + 1]);
    }

    // 3. Prevertex angles for circle 0.
    int b = 3 * m - 3;
    {
        const std::size_t K0 = m_circles[0].prevertexAngles.size();
        std::vector<double> psi(K0 - 1);
        for (std::size_t k = 0; k < psi.size(); ++k)
        {
            psi[k] = Xu[b + static_cast<int>(k)];
        }
        m_circles[0].prevertexAngles = mcsc::anglesFromUnconstrained(m_circles[0].prevertexAngles.front(), psi);
        b += static_cast<int>(psi.size());
    }

    // 4. Prevertex angles for circles 1..m-1.
    for (int j = 1; j < m; ++j)
    {
        const double theta1 = Xu[b];
        const std::size_t Kj = m_circles[j].prevertexAngles.size();
        std::vector<double> psi(Kj - 1);
        for (std::size_t k = 0; k < psi.size(); ++k)
        {
            psi[k] = Xu[b + 1 + static_cast<int>(k)];
        }
        m_circles[j].prevertexAngles = mcsc::anglesFromUnconstrained(theta1, psi);
        b += 1 + static_cast<int>(psi.size());
    }
}

bool MCSCCircleDomain::contains(const Complex& z) const
{
    // Inside the outer circle (component 0)...
    const auto& outer = m_circles[0];
    if (std::abs(z - outer.center) >= outer.radius)
    {
        return false;
    }

    // ...and outside every other (hole) circle.
    for (int j = 1; j < circleCount(); ++j)
    {
        const auto& hole = m_circles[j];
        if (std::abs(z - hole.center) <= hole.radius)
        {
            return false;
        }
    }

    return true;
}

void MCSCCircleDomain::transformBoundary(std::function<Complex(const Complex&)> transform)
{
    for (auto& circle : m_circles)
    {
        const Complex newCenter = transform(circle.center);
        const Complex probe = transform(circle.center + circle.radius);
        circle.center = newCenter;
        circle.radius = std::abs(probe - newCenter);
    }
}
