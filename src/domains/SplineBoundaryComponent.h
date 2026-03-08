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

#pragma once

#include "BoundaryComponent.h"

#include <vector>

/**
 * @brief Boundary component defined by a periodic cubic spline through control points.
 *
 * Ports the MATLAB spline_/interp_1/interp_2 algorithm (Hoskins & King 1972).
 * Parameterized by chordal arc-length, so totalLength() != 2*pi in general.
 *
 * Control points must form a closed curve (first == last, or closure is enforced).
 * Optional refinement upsamples the control points before computing the final spline.
 */
class SplineBoundaryComponent : public BoundaryComponent
{
public:
    /**
     * @brief Construct a spline boundary from control points
     * @param xpts X coordinates of control points (first == last for closure)
     * @param ypts Y coordinates of control points (first == last for closure)
     * @param refinement_N If > 0, upsample control points using interp_2 logic with
     *        nn = floor(refinement_N / (npts-1)) extra knots per interval
     */
    SplineBoundaryComponent(const std::vector<double>& xpts,
                            const std::vector<double>& ypts,
                            int refinement_N = 0);

    Complex evaluate(double s) const override;
    Complex evaluateDerivative(double s) const override;
    std::vector<Complex> sample(size_t numPoints) const override;
    double findParameterization(const Complex& z) const override;
    double totalLength() const override;

private:
    struct SplineCoefficients
    {
        std::vector<double> px, py;       // Control points (with closure)
        std::vector<double> x1, x2, x3;   // x spline derivative coefficients
        std::vector<double> y1, y2, y3;   // y spline derivative coefficients
        std::vector<double> h;            // Chordal segment distances
        double tl;                        // Total arc length
    };

    // Ported from spline_.m (Hoskins & King 1972)
    static SplineCoefficients computePeriodicSplineCoefficients(
        const std::vector<double>& x, const std::vector<double>& y);

    // Ported from interp_2.m — upsample control points
    static std::pair<std::vector<double>, std::vector<double>>
    upsampleControlPoints(const std::vector<double>& x,
                          const std::vector<double>& y, int nn);

    // Binary search for segment containing arc-length s
    int findSegment(double s) const;

    // Solve periodic tridiagonal system using Sherman-Morrison
    static std::vector<double> solvePeriodicTridiagonal(
        const std::vector<double>& a_sub,   // sub-diagonal (size n)
        const std::vector<double>& a_diag,  // diagonal (size n)
        const std::vector<double>& a_super, // super-diagonal (size n)
        double corner_bl,                    // bottom-left corner
        double corner_tr,                    // top-right corner
        const std::vector<double>& rhs);     // right-hand side (size n)

    // Solve tridiagonal system (Thomas algorithm)
    static std::vector<double> solveTridiagonal(
        std::vector<double> a_sub,
        std::vector<double> a_diag,
        std::vector<double> a_super,
        std::vector<double> rhs);

    std::vector<double> m_px, m_py;
    std::vector<double> m_x1, m_x2, m_x3;
    std::vector<double> m_y1, m_y2, m_y3;
    std::vector<double> m_h;
    std::vector<double> m_cumh;  // Cumulative h for segment lookup: [0, h[0], h[0]+h[1], ...]
    double m_tot_len;
    int m_num_segments;
};
