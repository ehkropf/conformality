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

#include "MCSCUnbounded.h"

#include "MCSCUnboundedObjectiveFunction.h"
#include "../core/ConformalMap.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
/// Wrap right into (left, left + 2*pi], matching extobjfun.m's/extmap.m's forward-arc convention.
double wrapForward(double left, double right)
{
    return (left > right) ? right + TWO_PI : right;
}
}  // namespace

MCSCUnbounded::MCSCUnbounded(int N, int ngj, MCSCContinuationSolver::Options continuationOptions)
    : m_N{N}
    , m_ngj{ngj}
    , m_continuationOptions{continuationOptions}
{
}

void MCSCUnbounded::compute(ConformalMap& map_instance, double /*target_accuracy*/)
{
    validateDomains(map_instance);

    auto targetDomain = std::dynamic_pointer_cast<MCSCPolygonalDomain>(map_instance.getTargetDomain());
    auto sourceDomain = std::dynamic_pointer_cast<MCSCCircleDomain>(map_instance.getSourceDomain());

    m_polygon.emplace(*targetDomain);

    MCSCUnboundedObjectiveFunction objective(*m_polygon, *sourceDomain, m_N, m_ngj);
    m_circle.emplace(objective.solve(m_continuationOptions));

    // objective.solve() doesn't surface the continuation solver's own residual, so recompute it
    // here at the solved point (cheap: one more evaluate() call) to report a real achieved
    // accuracy rather than a misleading hardcoded value.
    const Eigen::VectorXd finalResidual = objective.evaluate(m_circle->toUnconstrained());
    m_achieved_accuracy = finalResidual.lpNorm<Eigen::Infinity>();
    m_iteration_count = 0;  // MCSCContinuationSolver doesn't report a Newton-style iteration count.

    m_integrand.emplace(*m_polygon, *m_circle, m_N);

    std::vector<double> betaValues;
    for (int j = 0; j < m_polygon->getConnectivity(); ++j)
    {
        for (double a : m_polygon->getAlpha(j))
        {
            betaValues.push_back(1.0 - a);
        }
    }
    m_quadrature.emplace(betaValues, m_ngj);

    const int m = m_circle->circleCount();
    m_componentOffset.assign(m, 0);
    m_singularities.clear();
    for (int j = 0; j < m; ++j)
    {
        m_componentOffset[j] = static_cast<int>(m_singularities.size());
        const auto prevertices = m_circle->getPrevertices(j);
        m_singularities.insert(m_singularities.end(), prevertices.begin(), prevertices.end());
    }

    auto fprime = [this](const Complex& z) { return m_integrand->evalFPrime(z); };
    const auto& t0 = m_circle->getPrevertexAngles(0);
    const Complex Q12 = m_quadrature->integrateArc(
        t0[0], 0, wrapForward(t0[0], t0[1]), 1, m_circle->getCenter(0), m_circle->getRadius(0), fprime,
        m_singularities);
    if (!std::isfinite(std::abs(Q12)) || std::abs(Q12) < std::numeric_limits<double>::epsilon())
    {
        throw std::runtime_error(
            "MCSCUnbounded: normalization integral Q12 is degenerate (near zero or non-finite) -- "
            "target polygon vertices 0 and 1 on component 0 may be too close, or the parameter-problem "
            "solve did not converge to a valid circle configuration");
    }

    const auto& w0 = m_polygon->getVertices(0);
    m_A = (w0[1] - w0[0]) / Q12;
}

void MCSCUnbounded::ensureComputed() const
{
    if (!m_circle || !m_integrand || !m_quadrature || !m_polygon)
    {
        throw std::runtime_error("MCSCUnbounded: Map not computed yet");
    }
}

std::pair<int, double> MCSCUnbounded::nearestCircle(const Complex& z) const
{
    const int m = m_circle->circleCount();
    int best = 0;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (int j = 0; j < m; ++j)
    {
        const double d = std::abs(std::abs(z - m_circle->getCenter(j)) - m_circle->getRadius(j));
        if (d < bestDistance)
        {
            bestDistance = d;
            best = j;
        }
    }
    const double t = std::fmod(std::arg(z - m_circle->getCenter(best)) + TWO_PI, TWO_PI);
    return {best, t};
}

MCSCUnbounded::IntegrationPath MCSCUnbounded::integPath(const Complex& z) const
{
    const auto [j, t] = nearestCircle(z);
    const auto& tkj = m_circle->getPrevertexAngles(j);
    const int Kj = static_cast<int>(tkj.size());

    int k = 0;
    double minAbsDt = std::numeric_limits<double>::infinity();
    double minDt = 0.0;
    for (int i = 0; i < Kj; ++i)
    {
        double dt = tkj[i] - t;
        dt = std::fmod(dt + M_PI, TWO_PI);
        if (dt < 0.0)
        {
            dt += TWO_PI;
        }
        dt -= M_PI;
        if (std::abs(dt) < minAbsDt)
        {
            minAbsDt = std::abs(dt);
            minDt = dt;
            k = i;
        }
    }

    IntegrationPath path;
    path.circleIndex = j;
    path.nearestVertexFlat = m_componentOffset[j] + k;

    if (minAbsDt < 100.0 * std::numeric_limits<double>::epsilon())
    {
        path.arc = std::nullopt;
    }
    else
    {
        ArcSegment arc;
        arc.left = tkj[k];
        arc.leftVertex = path.nearestVertexFlat;
        arc.right = tkj[k] - minDt;
        arc.rightVertex = -1;
        path.arc = arc;
    }

    const double r = m_circle->getRadius(j);
    const double eps100 = 100.0 * std::numeric_limits<double>::epsilon() * std::max(r, 1.0);
    if (std::abs(z - m_circle->getCenter(j)) < r + eps100)
    {
        path.line = std::nullopt;
    }
    else
    {
        LineSegment line;
        line.left = m_circle->getCenter(j) + r * std::exp(Complex(0.0, 1.0) * t);
        line.leftVertex = path.arc ? -1 : path.nearestVertexFlat;
        path.line = line;
    }

    return path;
}

Complex MCSCUnbounded::map(const Complex& z) const
{
    ensureComputed();

    const auto path = integPath(z);
    const int j = path.circleIndex;

    // Decode (k, j) from the flat nearest-vertex index for the vertex-image lookup.
    int k = path.nearestVertexFlat - m_componentOffset[j];

    const auto& w = m_polygon->getVertices(j);
    Complex result = w[k];

    auto fprime = [this](const Complex& p) { return m_integrand->evalFPrime(p); };

    if (path.arc)
    {
        result += m_A
            * m_quadrature->integrateArc(
                path.arc->left, path.arc->leftVertex, path.arc->right, path.arc->rightVertex, m_circle->getCenter(j),
                m_circle->getRadius(j), fprime, m_singularities);
    }
    if (path.line)
    {
        result += m_A
            * m_quadrature->integrateLine(path.line->left, path.line->leftVertex, z, -1, fprime, m_singularities);
    }

    return result;
}

Complex MCSCUnbounded::inverseMap(const Complex& /*w*/) const
{
    throw std::runtime_error("MCSCUnbounded: inverseMap not implemented");
}

const MCSCCircleDomain& MCSCUnbounded::solvedCircleDomain() const
{
    ensureComputed();
    return *m_circle;
}

void MCSCUnbounded::validateSourceDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("MCSCUnbounded: Source domain cannot be null");
    }

    auto circle = std::dynamic_pointer_cast<MCSCCircleDomain>(domain);
    if (!circle)
    {
        throw std::invalid_argument("MCSCUnbounded: Source domain must be an MCSCCircleDomain");
    }
}

void MCSCUnbounded::validateTargetDomain(std::shared_ptr<Domain> domain) const
{
    if (!domain)
    {
        throw std::invalid_argument("MCSCUnbounded: Target domain cannot be null");
    }

    auto polygon = std::dynamic_pointer_cast<MCSCPolygonalDomain>(domain);
    if (!polygon)
    {
        throw std::invalid_argument("MCSCUnbounded: Target domain must be an MCSCPolygonalDomain");
    }

    if (!polygon->isUnbounded())
    {
        throw std::invalid_argument("MCSCUnbounded: Target domain must be unbounded");
    }
}
