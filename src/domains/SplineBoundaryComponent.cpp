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

#include "SplineBoundaryComponent.h"
#include "../numerics/RootFinder.h"
#include "../core/Tolerances.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

SplineBoundaryComponent::SplineBoundaryComponent(
    const std::vector<double>& xpts,
    const std::vector<double>& ypts,
    int refinement_N)
{
    if (xpts.size() != ypts.size())
    {
        throw std::invalid_argument("SplineBoundaryComponent: xpts and ypts must have the same size");
    }
    if (xpts.size() < 3)
    {
        throw std::invalid_argument("SplineBoundaryComponent: need at least 3 control points");
    }

    std::vector<double> px = xpts;
    std::vector<double> py = ypts;

    if (refinement_N > 0)
    {
        int npts = static_cast<int>(px.size());
        // Ensure closure before counting segments
        bool is_closed = (std::abs(px.front() - px.back()) <= SPLINE_CLOSURE_EPS
                       && std::abs(py.front() - py.back()) <= SPLINE_CLOSURE_EPS);
        int num_segments = is_closed ? npts - 1 : npts;
        int nn = refinement_N / num_segments;
        if (nn > 0)
        {
            auto [rx, ry] = upsampleControlPoints(px, py, nn);
            px = std::move(rx);
            py = std::move(ry);
        }
    }

    auto coeffs = computePeriodicSplineCoefficients(px, py);
    m_px = std::move(coeffs.px);
    m_py = std::move(coeffs.py);
    m_x1 = std::move(coeffs.x1);
    m_x2 = std::move(coeffs.x2);
    m_x3 = std::move(coeffs.x3);
    m_y1 = std::move(coeffs.y1);
    m_y2 = std::move(coeffs.y2);
    m_y3 = std::move(coeffs.y3);
    m_h = std::move(coeffs.h);
    m_tot_len = coeffs.tl;
    m_num_segments = static_cast<int>(m_px.size()) - 1;

    // Build cumulative arc-length breakpoints: [0, h[0], h[0]+h[1], ...]
    m_cumh.resize(m_num_segments + 1);
    m_cumh[0] = 0.0;
    for (int i = 0; i < m_num_segments; ++i)
    {
        m_cumh[i + 1] = m_cumh[i] + m_h[i];
    }
}

double SplineBoundaryComponent::totalLength() const
{
    return m_tot_len;
}

Complex SplineBoundaryComponent::evaluate(double s) const
{
    // Wrap s to [0, tot_len)
    s = s - std::floor(s / m_tot_len) * m_tot_len;

    int seg = findSegment(s);
    double ds = s - m_cumh[seg];

    // Cubic polynomial: p(ds) = p[i] + p1[i]*ds + p2[i]*ds^2/2 + p3[i]*ds^3/6
    // (matches mkpp coefficients in interp_1.m)
    double xval = m_px[seg] + m_x1[seg] * ds + m_x2[seg] * ds * ds / 2.0 + m_x3[seg] * ds * ds * ds / 6.0;
    double yval = m_py[seg] + m_y1[seg] * ds + m_y2[seg] * ds * ds / 2.0 + m_y3[seg] * ds * ds * ds / 6.0;

    return Complex(xval, yval);
}

Complex SplineBoundaryComponent::evaluateDerivative(double s) const
{
    s = s - std::floor(s / m_tot_len) * m_tot_len;

    int seg = findSegment(s);
    double ds = s - m_cumh[seg];

    // Derivative: p'(ds) = p1[i] + p2[i]*ds + p3[i]*ds^2/2
    double dxval = m_x1[seg] + m_x2[seg] * ds + m_x3[seg] * ds * ds / 2.0;
    double dyval = m_y1[seg] + m_y2[seg] * ds + m_y3[seg] * ds * ds / 2.0;

    return Complex(dxval, dyval);
}

std::vector<Complex> SplineBoundaryComponent::sample(size_t numPoints) const
{
    std::vector<Complex> samples;
    samples.reserve(numPoints);
    for (size_t j = 0; j < numPoints; ++j)
    {
        double s = m_tot_len * static_cast<double>(j) / static_cast<double>(numPoints);
        samples.push_back(evaluate(s));
    }
    return samples;
}

double SplineBoundaryComponent::findParameterization(const Complex& z) const
{
    auto objective = [this, &z](double s) -> double
    {
        Complex diff = this->evaluate(s) - z;
        return std::norm(diff);
    };

    try
    {
        double result = RootFinder::ternarySearch(objective, 0.0, m_tot_len, 1e-12);
        // Normalize to [0, tot_len)
        while (result < 0.0) result += m_tot_len;
        while (result >= m_tot_len) result -= m_tot_len;
        return result;
    }
    catch (const RootFinder::ConvergenceError&)
    {
        p_statusManager->reportWarning("SplineBoundaryComponent",
            "Ternary search failed to converge in findParameterization",
            "Falling back to coarse grid search");

        // Fallback: coarse grid search
        int n_search = 1000;
        double best_s = 0.0;
        double best_dist = std::numeric_limits<double>::max();
        for (int i = 0; i < n_search; ++i)
        {
            double s = m_tot_len * i / n_search;
            double dist = objective(s);
            if (dist < best_dist)
            {
                best_dist = dist;
                best_s = s;
            }
        }
        return best_s;
    }
}

int SplineBoundaryComponent::findSegment(double s) const
{
    // Binary search in cumh for the segment containing s
    // s is in [cumh[seg], cumh[seg+1])
    auto it = std::upper_bound(m_cumh.begin(), m_cumh.end(), s);
    int seg = static_cast<int>(it - m_cumh.begin()) - 1;
    if (seg < 0) seg = 0;
    if (seg >= m_num_segments) seg = m_num_segments - 1;
    return seg;
}

// ===== spline_.m port =====
SplineBoundaryComponent::SplineCoefficients
SplineBoundaryComponent::computePeriodicSplineCoefficients(
    const std::vector<double>& x_in, const std::vector<double>& y_in)
{
    std::vector<double> x = x_in;
    std::vector<double> y = y_in;

    // Ensure closure
    if (std::abs(x.front() - x.back()) > SPLINE_CLOSURE_EPS
        || std::abs(y.front() - y.back()) > SPLINE_CLOSURE_EPS)
    {
        x.push_back(x.front());
        y.push_back(y.front());
    }

    int n1 = static_cast<int>(x.size());
    int n = n1 - 1;

    // dx, dy, h
    std::vector<double> dx(n), dy(n), h(n);
    for (int i = 0; i < n; ++i)
    {
        dx[i] = x[i + 1] - x[i];
        dy[i] = y[i + 1] - y[i];
        h[i] = std::sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
    }

    // Validate no zero-length segments (prevents division by zero in coefficient computation)
    for (int i = 0; i < n; ++i)
    {
        if (h[i] < PIVOT_EPS)
        {
            throw std::invalid_argument(
                "SplineBoundaryComponent: zero-length chord at segment " + std::to_string(i)
                + " — duplicate consecutive control points detected");
        }
    }

    double tl = std::accumulate(h.begin(), h.end(), 0.0);

    // h(n1) = h(1) in MATLAB (0-indexed: h_ext[n] = h[0])
    // p = h(1:n), q = h(2:n1) = [h[1], ..., h[n-1], h[0]]
    std::vector<double> p(n), q(n);
    for (int i = 0; i < n; ++i)
    {
        p[i] = h[i];
        q[i] = h[(i + 1) % n];
    }

    // a = q./(p+q), b = 1-a
    std::vector<double> a(n), b(n);
    for (int i = 0; i < n; ++i)
    {
        a[i] = q[i] / (p[i] + q[i]);
        b[i] = 1.0 - a[i];
    }

    // Build the RHS vector d1 for x-derivatives:
    // d1 = 3*(a.*dx./p + b.*[dx(2:n); x(2)-x(n1)]./q)
    // In 0-indexed: dx_shifted[i] = dx[i+1] for i<n-1, dx_shifted[n-1] = x[1]-x[n1-1] = dx[0] (periodic)
    std::vector<double> d1_x(n), d1_y(n);
    for (int i = 0; i < n; ++i)
    {
        double dx_shifted = (i < n - 1) ? dx[i + 1] : (x[1] - x[n1 - 1]);
        double dy_shifted = (i < n - 1) ? dy[i + 1] : (y[1] - y[n1 - 1]);
        d1_x[i] = 3.0 * (a[i] * dx[i] / p[i] + b[i] * dx_shifted / q[i]);
        d1_y[i] = 3.0 * (a[i] * dy[i] / p[i] + b[i] * dy_shifted / q[i]);
    }

    // Build periodic tridiagonal system:
    // Diagonal: 2 (all entries)
    // Sub-diagonal: a[1], a[2], ..., a[n-1]  (MATLAB indexing: a(2:n) in positions -1)
    // Super-diagonal: b[0], b[1], ..., b[n-2] (MATLAB indexing: b(1:n-1) in positions +1)
    // Corner bottom-left: b[n-1]  (MATLAB: b(n) at position (-n+1))
    // Corner top-right: a[0]      (MATLAB: a(1) at position (n-1))
    std::vector<double> diag(n, 2.0);
    std::vector<double> sub(n), super(n);
    for (int i = 0; i < n; ++i)
    {
        sub[i] = (i > 0) ? a[i] : 0.0;     // sub[0] unused in Thomas, corner handles it
        super[i] = (i < n - 1) ? b[i] : 0.0;
    }
    double corner_bl = b[n - 1];   // ones(n-1,1);b(n) → last element of b
    double corner_tr = a[0];       // a(1);ones(n-1,1) → first element of a

    // Solve for x1 (first derivatives of x)
    auto x1_interior = solvePeriodicTridiagonal(sub, diag, super, corner_bl, corner_tr, d1_x);
    auto y1_interior = solvePeriodicTridiagonal(sub, diag, super, corner_bl, corner_tr, d1_y);

    // MATLAB: x1(2:n1) = x1; x1(1) = x1(n1);
    // So x1 has n1 elements with periodic wrap
    std::vector<double> x1(n1), y1(n1);
    for (int i = 0; i < n; ++i)
    {
        x1[i + 1] = x1_interior[i];
        y1[i + 1] = y1_interior[i];
    }
    x1[0] = x1[n];  // x1(1) = x1(n1) in MATLAB
    y1[0] = y1[n];

    // x2(2:n1) = 2*(x1(1:n) + 2*x1(2:n1) - 3*dx./p)./p
    std::vector<double> x2(n1), y2(n1);
    for (int i = 0; i < n; ++i)
    {
        x2[i + 1] = 2.0 * (x1[i] + 2.0 * x1[i + 1] - 3.0 * dx[i] / p[i]) / p[i];
        y2[i + 1] = 2.0 * (y1[i] + 2.0 * y1[i + 1] - 3.0 * dy[i] / p[i]) / p[i];
    }
    x2[0] = x2[n];
    y2[0] = y2[n];

    // x3 = diff(x2)./p, x3(n1) = x3(1)
    std::vector<double> x3(n1), y3(n1);
    for (int i = 0; i < n; ++i)
    {
        x3[i] = (x2[i + 1] - x2[i]) / p[i];
        y3[i] = (y2[i + 1] - y2[i]) / p[i];
    }
    x3[n] = x3[0];
    y3[n] = y3[0];

    // h(n1) = h(1) in MATLAB
    std::vector<double> h_ext(n1);
    for (int i = 0; i < n; ++i) h_ext[i] = h[i];
    h_ext[n] = h[0];

    return SplineCoefficients{x, y, x1, x2, x3, y1, y2, y3, h_ext, tl};
}

// ===== interp_2.m port =====
std::pair<std::vector<double>, std::vector<double>>
SplineBoundaryComponent::upsampleControlPoints(
    const std::vector<double>& x_in, const std::vector<double>& y_in, int nn)
{
    if (nn == 0)
    {
        return {x_in, y_in};
    }

    std::vector<double> x = x_in;
    std::vector<double> y = y_in;

    // Ensure closure
    if (std::abs(x.front() - x.back()) > SPLINE_CLOSURE_EPS
        || std::abs(y.front() - y.back()) > SPLINE_CLOSURE_EPS)
    {
        x.push_back(x.front());
        y.push_back(y.front());
    }

    int ns = static_cast<int>(x.size()) - 1;

    // Compute spline on coarse points
    auto coeffs = computePeriodicSplineCoefficients(x, y);

    // Build evaluation arc-lengths (interp_2.m logic):
    // mat1 = repmat((0:n-1)',1,ns)/n    → fractions [0, 1/nn, ..., (nn-1)/nn]
    // mat2 = repmat(h(1:ns)',n,1)        → segment lengths
    // v = cumsum([0;h(1:ns-1)])           → cumulative segment starts
    // s = mat3 + mat1.*mat2               → fine arc-lengths
    std::vector<double> s;
    s.reserve(ns * nn + 1);
    double cum = 0.0;
    for (int seg = 0; seg < ns; ++seg)
    {
        for (int k = 0; k < nn; ++k)
        {
            double frac = static_cast<double>(k) / static_cast<double>(nn);
            s.push_back(cum + frac * coeffs.h[seg]);
        }
        cum += coeffs.h[seg];
    }
    // s(ns*n+1) = 0 in MATLAB — wraps back to start
    s.push_back(0.0);

    // Build cumulative breakpoints for ppval
    int n_segs = static_cast<int>(coeffs.px.size()) - 1;
    std::vector<double> cumh(n_segs + 1);
    cumh[0] = 0.0;
    for (int i = 0; i < n_segs; ++i)
    {
        cumh[i + 1] = cumh[i] + coeffs.h[i];
    }

    // Evaluate spline at each s (interp_1.m logic)
    std::vector<double> fx, fy;
    fx.reserve(s.size());
    fy.reserve(s.size());

    for (double sv : s)
    {
        // Wrap to [0, tl)
        sv = sv - std::floor(sv / coeffs.tl) * coeffs.tl;

        // Find segment
        auto it = std::upper_bound(cumh.begin(), cumh.end(), sv);
        int seg = static_cast<int>(it - cumh.begin()) - 1;
        if (seg < 0) seg = 0;
        if (seg >= n_segs) seg = n_segs - 1;

        double ds = sv - cumh[seg];
        double xval = coeffs.px[seg] + coeffs.x1[seg] * ds + coeffs.x2[seg] * ds * ds / 2.0
                    + coeffs.x3[seg] * ds * ds * ds / 6.0;
        double yval = coeffs.py[seg] + coeffs.y1[seg] * ds + coeffs.y2[seg] * ds * ds / 2.0
                    + coeffs.y3[seg] * ds * ds * ds / 6.0;
        fx.push_back(xval);
        fy.push_back(yval);
    }

    return {fx, fy};
}

// ===== Tridiagonal solvers =====

std::vector<double> SplineBoundaryComponent::solveTridiagonal(
    std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
    int n = static_cast<int>(b.size());
    // Forward sweep
    for (int i = 1; i < n; ++i)
    {
        if (std::abs(b[i - 1]) < PIVOT_EPS)
        {
            throw std::runtime_error(
                "SplineBoundaryComponent: tridiagonal solver encountered zero pivot at index "
                + std::to_string(i - 1));
        }
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }
    // Back substitution
    std::vector<double> x(n);
    if (std::abs(b[n - 1]) < PIVOT_EPS)
    {
        throw std::runtime_error(
            "SplineBoundaryComponent: tridiagonal solver encountered zero pivot at index "
            + std::to_string(n - 1));
    }
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    return x;
}

std::vector<double> SplineBoundaryComponent::solvePeriodicTridiagonal(
    const std::vector<double>& a_sub,
    const std::vector<double>& a_diag,
    const std::vector<double>& a_super,
    double corner_bl,
    double corner_tr,
    const std::vector<double>& rhs)
{
    int n = static_cast<int>(a_diag.size());

    // Sherman-Morrison: A = A' + u*v^T
    // where A' is tridiagonal (with modified corners) and u*v^T accounts for periodic terms
    //
    // gamma = -a_diag[0] (arbitrary nonzero)
    double gamma = -a_diag[0];

    // Modified diagonal: a_diag[0] - gamma, a_diag[n-1] - corner_bl*corner_tr/gamma
    std::vector<double> diag_mod = a_diag;
    diag_mod[0] -= gamma;
    diag_mod[n - 1] -= corner_bl * corner_tr / gamma;

    // u = [gamma, 0, ..., 0, corner_bl]
    // v = [1, 0, ..., 0, corner_tr/gamma]
    std::vector<double> u(n, 0.0), v(n, 0.0);
    u[0] = gamma;
    u[n - 1] = corner_bl;
    v[0] = 1.0;
    v[n - 1] = corner_tr / gamma;

    // Solve A'*y = rhs and A'*z = u
    auto y = solveTridiagonal(
        std::vector<double>(a_sub.begin(), a_sub.end()),
        diag_mod,
        std::vector<double>(a_super.begin(), a_super.end()),
        std::vector<double>(rhs.begin(), rhs.end()));

    auto z = solveTridiagonal(
        std::vector<double>(a_sub.begin(), a_sub.end()),
        diag_mod,
        std::vector<double>(a_super.begin(), a_super.end()),
        u);

    // x = y - z * (v^T * y) / (1 + v^T * z)
    double vTy = 0.0, vTz = 0.0;
    for (int i = 0; i < n; ++i)
    {
        vTy += v[i] * y[i];
        vTz += v[i] * z[i];
    }

    double denom = 1.0 + vTz;
    if (std::abs(denom) < GEOMETRIC_COINCIDENCE_EPS)
    {
        throw std::runtime_error(
            "SplineBoundaryComponent: periodic tridiagonal system is singular "
            "(degenerate control point geometry)");
    }
    double factor = vTy / denom;
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i)
    {
        x[i] = y[i] - factor * z[i];
    }

    return x;
}
