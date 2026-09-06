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

#include <Eigen/Dense>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * @brief Gauss-Jacobi quadrature with vertex-singularity handling, for the MCSC reflection
 *        method (dissertation Sec 3.4; port of gjquad.m).
 *
 * Evaluates integrals of MCSC's reflection-method integrand (e.g. MCSCReflectionIntegrand's
 * f'(z), #163) along a path built from straight-line and circular-arc segments, where the
 * integrand has algebraic endpoint singularities at polygon-vertex preimages with strength
 * beta_{k,j} = 1 - alpha_{k,j}. Gauss-Jacobi quadrature handles the singular endpoint of a
 * segment directly (by baking the singularity's strength into the quadrature weight function);
 * the "one-half rule" (dissertation Sec 3.4: no singularity may lie closer to an integration
 * subinterval than half that subinterval's length) adaptively subdivides a segment so that
 * *other* nearby singularities (original or reflected prevertices) don't corrupt accuracy.
 *
 * Known limitation, faithfully ported from gjquad.m's arcq/lineq: the one-half-rule check is
 * evaluated incrementally from the current walking position, not as a lookahead over the whole
 * candidate subinterval -- it reliably subdivides away from a singularity near wherever the
 * walk currently stands (i.e. near a path endpoint, which is the case that matters for MCSC,
 * since nearby reflected prevertices cluster near path endpoints), but a singularity sitting
 * near the *interior* of an otherwise-long candidate step may not trigger any subdivision on
 * that step at all. This mirrors the MATLAB reference exactly (`abs(a - sng)` in both
 * `arcq`/`lineq`), not a defect introduced by this port.
 *
 * Deliberately decoupled from any specific integrand type: the integrand is passed as a
 * std::function callable and the singularity list as an explicit argument to each integrate
 * call (matching RootFinder's callable-based style), rather than mirroring gjquad.m's coupling
 * to a specific Q.fp integrand object with its own .sing property. This keeps quadrature
 * reusable across the eventual bounded-case integrand too, and independently testable against
 * known closed-form integrals without needing an MCSCReflectionIntegrand fixture.
 *
 * gjquad.m's quadrature (marker superclass, no behavior) and pvnum (pure (k,j) -> flat-index
 * arithmetic) are not ported: callers do that index-flattening themselves when building the
 * flat betaValues list passed to the constructor.
 */
class GaussJacobiQuadrature
{
public:
    /**
     * @brief Exception thrown when the one-half-rule subdivision loop exceeds its subdivision
     *        cap without reaching the far endpoint (port of gjquad.m's 'mcsc:manysubdiv' error).
     */
    class ConvergenceError : public std::runtime_error
    {
    public:
        explicit ConvergenceError(const std::string& message)
            : std::runtime_error(message)
        {
        }
    };

    /// Integrand callable: evaluates the function being integrated at a point.
    using Integrand = std::function<Complex(const Complex&)>;

    /**
     * @brief Compute Gauss-Jacobi nodes and weights on [-1,1] for the weight function
     *        (1-x)^alf * (1+x)^bet (port of gjquad.m's gaussj(), Golub-Welsch algorithm).
     *
     * Builds the symmetric tridiagonal "Ritz matrix" from the closed-form Jacobi three-term
     * Lanczos recurrence coefficients; its eigenvalues are the quadrature nodes and its
     * eigenvectors' first components (scaled) are the weights.
     *
     * @param n Number of nodes (n >= 1).
     * @param alf Jacobi weight exponent for (1-x) (alf > -1).
     * @param bet Jacobi weight exponent for (1+x) (bet > -1).
     * @return (nodes, weights), both size n, nodes sorted ascending.
     * @throws std::invalid_argument if n < 1, or alf <= -1, or bet <= -1.
     */
    static std::pair<Eigen::VectorXd, Eigen::VectorXd> gaussianJacobiNodes(int n, double alf, double bet);

    /**
     * @brief Build the per-vertex Gauss-Jacobi node/weight tables (port of gjquad.m's
     *        calc_qdata()).
     *
     * @param betaValues Flat per-vertex beta_{k,j} = 1 - alpha_{k,j} values (caller flattens
     *        (k,j) pairs into a single index, e.g. via MCSCPolygonalDomain's per-component
     *        alpha lists concatenated in component order). Vertices with beta <= -1 have no
     *        valid Gauss-Jacobi weight function and are skipped (not usable as a singular
     *        endpoint in integrateArc/integrateLine -- see their docs).
     * @param ngj Node count per column, floored to 4 (mirrors MATLAB's max(ngj,4)).
     * @param useHalfRule Whether to apply the one-half-rule adaptive subdivision (true, the
     *        default) or integrate each segment in a single Gauss-Jacobi/ordinary step
     *        (false) -- mirrors gjquad.m's halfrule flag.
     */
    explicit GaussJacobiQuadrature(std::vector<double> betaValues, int ngj = 12, bool useHalfRule = true);

    /**
     * @brief Integrate f along a circular arc, center c radius r, from angle left to angle
     *        right (port of gjquad.m's arcq()).
     *
     * @param left Starting angle (radians).
     * @param leftVertex Flat vertex index (as passed to the constructor's betaValues) if the
     *        left endpoint is a genuine prevertex singularity, or -1 if not.
     * @param right Ending angle (radians).
     * @param rightVertex Same convention as leftVertex, for the right endpoint. If rightVertex
     *        is a prevertex, the interval is split at its midpoint and integrated inward from
     *        both singular ends (matches arcq's rpn-triggered split).
     * @param c Arc center.
     * @param r Arc radius (r > 0).
     * @param f Integrand callable.
     * @param singularities Points (original or reflected prevertices) tracked by the one-half
     *        rule when subdividing -- does not need to include the endpoints themselves.
     * @return The integral's value.
     * @throws std::invalid_argument if r <= 0, or leftVertex/rightVertex references a vertex
     *         with no valid node/weight column (beta <= -1).
     * @throws ConvergenceError if more than 100 subdivisions are needed to satisfy the
     *         one-half rule for some sub-segment.
     */
    Complex integrateArc(
        double left,
        int leftVertex,
        double right,
        int rightVertex,
        const Complex& c,
        double r,
        const Integrand& f,
        const std::vector<Complex>& singularities
    ) const;

    /**
     * @brief Integrate f along a straight line segment from left to right (port of gjquad.m's
     *        lineq()). Same leftVertex/rightVertex/singularities conventions as integrateArc.
     * @throws std::invalid_argument if leftVertex/rightVertex references a vertex with no valid
     *         node/weight column (beta <= -1).
     * @throws ConvergenceError if more than 100 subdivisions are needed to satisfy the
     *         one-half rule for some sub-segment.
     */
    Complex integrateLine(
        const Complex& left,
        int leftVertex,
        const Complex& right,
        int rightVertex,
        const Integrand& f,
        const std::vector<Complex>& singularities
    ) const;

private:
    int m_ngj;
    bool m_useHalfRule;
    std::vector<double> m_betaValues;                // as passed to the constructor
    std::vector<bool> m_hasVertexTable;              // whether betaValues[k] > -1 (has a valid column)
    std::vector<Eigen::VectorXd> m_vertexNodes;       // per-vertex node columns, empty if !m_hasVertexTable[k]
    std::vector<Eigen::VectorXd> m_vertexWeights;     // per-vertex weight columns
    Eigen::VectorXd m_ordinaryNodes;                  // beta=0 column, always present
    Eigen::VectorXd m_ordinaryWeights;

    void validateVertexIndex(int vertex) const;

    /**
     * @brief Shared one-half-rule step-halving/subdivision loop used by both integrateArc and
     *        integrateLine.
     *
     * The one-half rule's step-length bookkeeping needs true Euclidean (arc-length) distances,
     * not raw parameter differences -- for integrateArc, a parameter difference of dtheta
     * corresponds to an arc length of r*dtheta, so arcLengthScale = r converts between them
     * (integrateLine's parameter *is* already the point itself, so arcLengthScale = 1).
     *
     * @param leftParam Parameter value (angle for arc, the point itself -- as a real multiple
     *        of (right-left) -- is not used; see pointAt) at the left endpoint.
     * @param rightParam Parameter value at the right endpoint.
     * @param leftVertex/rightVertex As in integrateArc/integrateLine.
     * @param arcLengthScale Converts a parameter-space distance into a true Euclidean distance
     *        (r for arcs, 1 for lines).
     * @param pointAt Maps a parameter value to the corresponding point in the complex plane.
     * @param jacobianAt Maps a parameter value to d(point)/d(parameter) at that value (arc:
     *        i*r*exp(i*theta); line: (right - left), constant).
     * @param f Integrand callable.
     * @param singularities Tracked singularity points for the one-half rule.
     * @throws ConvergenceError if more than 100 subdivisions are needed.
     */
    Complex integrateWithHalfRule(
        double leftParam,
        int leftVertex,
        double rightParam,
        int rightVertex,
        double arcLengthScale,
        const std::function<Complex(double)>& pointAt,
        const std::function<Complex(double)>& jacobianAt,
        const Integrand& f,
        const std::vector<Complex>& singularities
    ) const;
};
