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

#include "MCSCUnboundedObjectiveFunction.h"

#include "../core/Types.h"

#include <cmath>
#include <stdexcept>

namespace
{
/// Wrap right into (left, left + 2*pi] so an arc integral always runs forward, matching
/// extobjfun.m's `right(left>right) = right(left>right) + 2*pi;`.
double wrapForward(double left, double right)
{
    return (left > right) ? right + TWO_PI : right;
}
}

GaussJacobiQuadrature MCSCUnboundedObjectiveFunction::makeQuadrature(const MCSCPolygonalDomain& polygon, int ngj)
{
    std::vector<double> betaValues;
    for (int j = 0; j < polygon.getConnectivity(); ++j)
    {
        const auto& alpha = polygon.getAlpha(j);
        for (double a : alpha)
        {
            betaValues.push_back(1.0 - a);
        }
    }
    return GaussJacobiQuadrature(std::move(betaValues), ngj);
}

MCSCUnboundedObjectiveFunction::MCSCUnboundedObjectiveFunction(
    const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& initialCircle, int N, int ngj
)
    : m_polygon{polygon}
    , m_N{N}
    , m_circle{initialCircle}
    , m_integrand{polygon, m_circle, N}
    , m_quadrature{makeQuadrature(polygon, ngj)}
{
}

Eigen::VectorXd MCSCUnboundedObjectiveFunction::evaluate(const Eigen::VectorXd& Xu)
{
    m_circle.setFromUnconstrained(Xu);
    m_integrand.rebuild(m_polygon, m_circle, m_N);

    const int m = m_circle.circleCount();

    // Flat (0-based) vertex index of vertex k on component j, matching pvnum.m's convention.
    std::vector<int> componentOffset(m, 0);
    for (int j = 1; j < m; ++j)
    {
        componentOffset[j] = componentOffset[j - 1] + static_cast<int>(m_circle.prevertexCount(j - 1));
    }

    auto fprime = [this](const Complex& z) { return m_integrand.evalFPrime(z); };

    // Scaling/rotation constant A (extobjfun.m: Q12 = arcq(...); A = (vl(2,1)-vl(1,1))/Q12).
    const auto& t0 = m_circle.getPrevertexAngles(0);
    const Complex Q12 = m_quadrature.integrateArc(
        t0[0], 0, wrapForward(t0[0], t0[1]), 1, m_circle.getCenter(0), m_circle.getRadius(0), fprime, {});
    const auto& w0 = m_polygon.getVertices(0);
    const Complex A = (w0[1] - w0[0]) / Q12;

    // Side lengths SL(k, j): arc integral around each consecutive prevertex pair on circle j.
    std::vector<std::vector<Complex>> SL(m);
    for (int j = 0; j < m; ++j)
    {
        const auto& t = m_circle.getPrevertexAngles(j);
        const int Kj = static_cast<int>(t.size());
        SL[j].resize(Kj);
        for (int k = 0; k < Kj; ++k)
        {
            const int kNext = (k + 1) % Kj;
            const double left = t[k];
            const double rightAngle = wrapForward(left, t[kNext]);
            const int leftFlat = componentOffset[j] + k;
            const int rightFlat = componentOffset[j] + kNext;
            SL[j][k] = A
                * m_quadrature.integrateArc(
                    left, leftFlat, rightAngle, rightFlat, m_circle.getCenter(j), m_circle.getRadius(j), fprime, {});
        }
    }

    // Position T(j-1): line integral from z(0,0) to z(0,j), for j = 1, ..., m-1.
    const Complex z00 = m_circle.getPrevertices(0)[0];
    std::vector<Complex> T(m > 0 ? m - 1 : 0);
    for (int j = 1; j < m; ++j)
    {
        const Complex z0j = m_circle.getPrevertices(j)[0];
        // Right-endpoint flat vertex index: vertex 0 of component j (extobjfun.m computes this
        // same value via `1 + cumsum(vc(1:j-1))'` with j left at m after the preceding
        // side-length loop -- that cumsum is a full (m-1)-vector, so despite the stale scalar j
        // it still lands on the correct per-position index; componentOffset[j] is the direct
        // 0-based equivalent).
        T[j - 1] = A * m_quadrature.integrateLine(z00, 0, z0j, componentOffset[j], fprime, {});
    }

    Eigen::VectorXd F(Xu.size());

    // 1. Orientation + first-side-length conditions for circles 1..m-1.
    for (int j = 1; j < m; ++j)
    {
        const Complex diff = SL[j][0] - (m_polygon.getVertices(j)[1] - m_polygon.getVertices(j)[0]);
        F[2 * (j - 1)] = diff.real();
        F[2 * (j - 1) + 1] = diff.imag();
    }

    // 2. Position conditions for circles 1..m-1 relative to circle 0.
    const int posBase = 2 * (m - 1);
    for (int j = 1; j < m; ++j)
    {
        const Complex diff = T[j - 1] - (m_polygon.getVertices(j)[0] - m_polygon.getVertices(0)[0]);
        F[posBase + 2 * (j - 1)] = diff.real();
        F[posBase + 2 * (j - 1) + 1] = diff.imag();
    }

    // 3. Remaining side-length conditions (magnitude only), for every circle.
    int idx = 4 * (m - 1);
    for (int j = 0; j < m; ++j)
    {
        const auto& w = m_polygon.getVertices(j);
        const int Kj = static_cast<int>(w.size());
        // Circle 0 skips k=0 (fixed by A); circles 1..m-1 skip k=0 (covered by step 1). Either
        // way the loop starts at k=1 (matches extobjfun.m's `for k = 2:vc(j)-1` plus its final
        // wraparound term, adjusted to 0-based).
        for (int k = 1; k < Kj - 1; ++k)
        {
            F[idx++] = std::abs(SL[j][k]) - std::abs(w[k + 1] - w[k]);
        }
        F[idx++] = std::abs(SL[j][Kj - 1]) - std::abs(w[0] - w[Kj - 1]);
    }

    // A degenerate circle configuration (e.g. Q12 collapsing toward zero during an exploratory
    // continuation step) can silently turn the scaling constant A -- and therefore every entry of
    // F -- into Inf/NaN with no exception raised along the way (GaussJacobiQuadrature's own
    // finiteness guard only covers its internal per-node arithmetic, not this function's
    // subsequent division by Q12). Surface that explicitly rather than returning a non-finite
    // residual that would otherwise defeat MCSCContinuationSolver's own tolerance checks.
    if (!F.allFinite())
    {
        throw std::runtime_error(
            "MCSCUnboundedObjectiveFunction::evaluate: non-finite residual (degenerate circle configuration)");
    }

    return F;
}

MCSCCircleDomain MCSCUnboundedObjectiveFunction::solve(MCSCContinuationSolver::Options options)
{
    const Eigen::VectorXd x0 = m_circle.toUnconstrained();

    MCSCContinuationSolver solver(
        [this](const Eigen::VectorXd& Xu) { return evaluate(Xu); }, options);
    const auto result = solver.solve(x0);

    m_circle.setFromUnconstrained(result.solution);
    return m_circle;
}
