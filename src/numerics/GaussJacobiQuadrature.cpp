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

#include "GaussJacobiQuadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
constexpr int MAX_SUBDIVISIONS = 100;

double distanceToNearestSingularity(const Complex& point, const std::vector<Complex>& singularities)
{
    double minDist = std::numeric_limits<double>::infinity();
    for (const Complex& s : singularities)
    {
        minDist = std::min(minDist, std::abs(point - s));
    }
    return minDist;
}
}  // namespace

std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussJacobiQuadrature::gaussianJacobiNodes(int n, double alf, double bet)
{
    if (n < 1)
    {
        throw std::invalid_argument("GaussJacobiQuadrature::gaussianJacobiNodes: n must be at least 1");
    }
    if (alf <= -1.0 || bet <= -1.0)
    {
        throw std::invalid_argument("GaussJacobiQuadrature::gaussianJacobiNodes: alf and bet must be > -1");
    }

    const double apb = alf + bet;

    // Three-term Lanczos recurrence coefficients for the Jacobi weight (closed form).
    Eigen::VectorXd a(n);
    Eigen::VectorXd b(n > 1 ? n - 1 : 0);

    a[0] = (bet - alf) / (apb + 2.0);
    if (n > 1)
    {
        b[0] = std::sqrt(4.0 * (1.0 + alf) * (1.0 + bet) / ((apb + 3.0) * (apb + 2.0) * (apb + 2.0)));
    }
    for (int k = 2; k <= n; ++k)
    {
        const double kd = static_cast<double>(k);
        a[k - 1] = apb * (bet - alf) / ((apb + 2.0 * kd) * (apb + 2.0 * kd - 2.0));
    }
    for (int k = 2; k <= n - 1; ++k)
    {
        const double kd = static_cast<double>(k);
        b[k - 1] = std::sqrt(
            4.0 * kd * (kd + alf) * (kd + bet) * (kd + apb) /
            (((apb + 2.0 * kd) * (apb + 2.0 * kd) - 1.0) * (apb + 2.0 * kd) * (apb + 2.0 * kd)));
    }

    Eigen::VectorXd nodes(n);
    Eigen::VectorXd weights(n);

    if (n > 1)
    {
        Eigen::MatrixXd ritz = Eigen::MatrixXd::Zero(n, n);
        for (int k = 0; k < n; ++k)
        {
            ritz(k, k) = a[k];
        }
        for (int k = 0; k < n - 1; ++k)
        {
            ritz(k, k + 1) = b[k];
            ritz(k + 1, k) = b[k];
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(ritz);
        nodes = solver.eigenvalues();

        const double c =
            std::pow(2.0, apb + 1.0) * std::tgamma(alf + 1.0) * std::tgamma(bet + 1.0) / std::tgamma(apb + 2.0);
        for (int k = 0; k < n; ++k)
        {
            const double v0 = solver.eigenvectors()(0, k);
            weights[k] = c * v0 * v0;
        }
    }
    else
    {
        nodes[0] = a[0];
        weights[0] =
            std::pow(2.0, apb + 1.0) * std::tgamma(alf + 1.0) * std::tgamma(bet + 1.0) / std::tgamma(apb + 2.0);
    }

    // SelfAdjointEigenSolver already returns eigenvalues in ascending order, matching gaussj.m's
    // explicit sort -- nothing further to do here.
    return {nodes, weights};
}

GaussJacobiQuadrature::GaussJacobiQuadrature(std::vector<double> betaValues, int ngj, bool useHalfRule)
    : m_ngj(std::max(ngj, 4))
    , m_useHalfRule(useHalfRule)
    , m_betaValues(std::move(betaValues))
{
    const auto [ordinaryNodes, ordinaryWeights] = gaussianJacobiNodes(m_ngj, 0.0, 0.0);
    m_ordinaryNodes = ordinaryNodes;
    m_ordinaryWeights = ordinaryWeights;

    m_hasVertexTable.resize(m_betaValues.size());
    m_vertexNodes.resize(m_betaValues.size());
    m_vertexWeights.resize(m_betaValues.size());
    for (std::size_t k = 0; k < m_betaValues.size(); ++k)
    {
        if (m_betaValues[k] <= -1.0)
        {
            m_hasVertexTable[k] = false;
            continue;
        }
        const auto [nodes, weights] = gaussianJacobiNodes(m_ngj, 0.0, m_betaValues[k]);
        m_vertexNodes[k] = nodes;
        m_vertexWeights[k] = weights;
        m_hasVertexTable[k] = true;
    }
}

void GaussJacobiQuadrature::validateVertexIndex(int vertex) const
{
    if (vertex < 0)
    {
        return;
    }
    if (vertex >= static_cast<int>(m_hasVertexTable.size()) || !m_hasVertexTable[static_cast<std::size_t>(vertex)])
    {
        throw std::invalid_argument(
            "GaussJacobiQuadrature: vertex index has no valid Gauss-Jacobi table (beta <= -1 or out of range)");
    }
}

Complex GaussJacobiQuadrature::integrateWithHalfRule(
    double leftParam,
    int leftVertex,
    double rightParam,
    int rightVertex,
    double arcLengthScale,
    const std::function<Complex(double)>& pointAt,
    const std::function<Complex(double)>& jacobianAt,
    const Integrand& f,
    const std::vector<Complex>& singularities
) const
{
    validateVertexIndex(leftVertex);
    validateVertexIndex(rightVertex);

    // If the right endpoint is itself a singularity, split symmetrically at the midpoint and
    // integrate both halves inward from their own singular end -- port of arcq/lineq's
    // rpn-triggered split. Integrating (right -> mid) and negating gives (mid -> right); adding
    // that to (left -> mid) gives the full (left -> right) integral.
    if (rightVertex >= 0)
    {
        const double mid = 0.5 * (leftParam + rightParam);
        Complex leftHalf = integrateWithHalfRule(
            leftParam, leftVertex, mid, -1, arcLengthScale, pointAt, jacobianAt, f, singularities);
        Complex rightHalf = integrateWithHalfRule(
            rightParam, rightVertex, mid, -1, arcLengthScale, pointAt, jacobianAt, f, singularities);
        return leftHalf - rightHalf;
    }

    const double dir = (rightParam > leftParam) ? 1.0 : -1.0;
    const double totalLength = arcLengthScale * std::abs(rightParam - leftParam);

    Complex accumulated{0.0, 0.0};
    double a = leftParam;
    int subdivisions = -1;
    bool firstStep = true;

    while (arcLengthScale * std::abs(leftParam - a) + 10.0 * std::numeric_limits<double>::epsilon() * totalLength <
           totalLength)
    {
        ++subdivisions;
        if (subdivisions > MAX_SUBDIVISIONS)
        {
            throw ConvergenceError(
                "GaussJacobiQuadrature: too many subdivisions detected while enforcing the one-half rule");
        }

        double b;
        if (m_useHalfRule)
        {
            const double distToFarEnd = arcLengthScale * std::abs(a - rightParam);
            const double distToSingularity = distanceToNearestSingularity(pointAt(a), singularities);
            const double len = std::min(distToFarEnd, 2.0 * distToSingularity);
            b = a + dir * (len / arcLengthScale);
        }
        else
        {
            b = rightParam;
        }

        const bool useVertexTable = firstStep && leftVertex >= 0;
        const Eigen::VectorXd& nodes = useVertexTable ? m_vertexNodes[static_cast<std::size_t>(leftVertex)]
                                                       : m_ordinaryNodes;
        const Eigen::VectorXd& weights = useVertexTable ? m_vertexWeights[static_cast<std::size_t>(leftVertex)]
                                                         : m_ordinaryWeights;
        const double beta = useVertexTable ? m_betaValues[static_cast<std::size_t>(leftVertex)] : 0.0;

        Complex stepValue{0.0, 0.0};
        for (int i = 0; i < nodes.size(); ++i)
        {
            const double x = nodes[i];
            const double w = weights[i];
            const double t = a + (b - a) * (x + 1.0) / 2.0;
            const Complex z = pointAt(t);
            const Complex fz = f(z);
            const Complex jac = jacobianAt(t);

            if (useVertexTable)
            {
                // Gauss-Jacobi weight already accounts for the (1+x)^beta singularity factor at
                // the left endpoint (x = -1); the integrand's own (z - z_vertex)^beta factor
                // must be divided out to avoid double-counting it (matches arcq/lineq's
                // division by (1+x)^betav(lpn)).
                stepValue += w * (0.5 * (b - a) * fz * jac) / std::pow(1.0 + x, beta);
            }
            else
            {
                stepValue += w * (0.5 * (b - a) * fz * jac);
            }
        }
        accumulated += stepValue;

        a = b;
        firstStep = false;
    }

    return accumulated;
}

Complex GaussJacobiQuadrature::integrateArc(
    double left,
    int leftVertex,
    double right,
    int rightVertex,
    const Complex& c,
    double r,
    const Integrand& f,
    const std::vector<Complex>& singularities
) const
{
    if (r <= 0.0)
    {
        throw std::invalid_argument("GaussJacobiQuadrature::integrateArc: radius must be positive");
    }

    auto pointAt = [c, r](double theta) { return c + r * std::exp(Complex(0.0, theta)); };
    auto jacobianAt = [r](double theta) { return Complex(0.0, r) * std::exp(Complex(0.0, theta)); };

    return integrateWithHalfRule(left, leftVertex, right, rightVertex, r, pointAt, jacobianAt, f, singularities);
}

Complex GaussJacobiQuadrature::integrateLine(
    const Complex& left,
    int leftVertex,
    const Complex& right,
    int rightVertex,
    const Integrand& f,
    const std::vector<Complex>& singularities
) const
{
    const Complex direction = right - left;

    // Parametrize by true arc length along the segment (t in [0, |right-left|]) so the shared
    // one-half-rule helper's length comparisons are already in Euclidean units (arcLengthScale
    // = 1), matching integrateArc's convention.
    const double length = std::abs(direction);
    auto pointAt = [left, direction, length](double t) {
        return length > 0.0 ? left + (t / length) * direction : left;
    };
    auto jacobianAt = [direction, length](double) { return length > 0.0 ? direction / length : Complex(0.0, 0.0); };

    return integrateWithHalfRule(0.0, leftVertex, length, rightVertex, 1.0, pointAt, jacobianAt, f, singularities);
}
