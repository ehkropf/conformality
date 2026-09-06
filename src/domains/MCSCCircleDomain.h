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

#include "../core/Types.h"
#include "Domain.h"

#include <Eigen/Dense>
#include <vector>

/**
 * @brief Canonical (source) circle domain for the MCSC method (port of circdomain.m).
 *
 * Represents an m-connected circle domain: m circles with centers c_j, radii r_j, and per-circle
 * prevertex angles t(1,j) < t(2,j) < ... < t(Kj,j). By convention (matching circdomain.m and the
 * dissertation), component 0 is the outer circle with c_0 = 0, r_0 = 1, and t(1,0) fixed at 0;
 * components 1, ..., m-1 are the remaining circles (interior circles for the bounded case, or the
 * other boundary circles for the unbounded case).
 *
 * This is a pure data representation plus the parameter-vector packing used by the MCSC Newton
 * parameter-problem solve (#165, #168) -- no mapping, quadrature, or reflection logic lives here.
 *
 * Deliberately a sibling of FornbergCanonicalDomain, not a subclass: FornbergCanonicalDomain
 * carries Boundary-component storage and an FFT-oriented sample count that are meaningless here,
 * since the reflection method enforces conditions at discrete prevertex points, not via dense
 * boundary sampling.
 */
class MCSCCircleDomain : public Domain
{
public:
    /**
     * @brief One circle's data: center, radius, and ordered prevertex angles.
     */
    struct CircleData
    {
        Complex center;
        double radius;
        std::vector<double> prevertexAngles;
    };

    /**
     * @brief Construct an m-connected circle domain from per-circle data.
     * @param circles Per-circle center/radius/prevertex-angle data (size m >= 1). Component 0 is
     *                treated as the outer circle by convention.
     * @throws std::invalid_argument if circles is empty, any radius is non-positive, any circle
     *         has fewer than 2 prevertices, or any circle's prevertex angles are not strictly
     *         increasing on [t(1,j), t(1,j) + 2*pi).
     */
    explicit MCSCCircleDomain(std::vector<CircleData> circles);

    ~MCSCCircleDomain() override = default;

    /**
     * @brief Number of circles (== getConnectivity()).
     */
    int circleCount() const
    {
        return static_cast<int>(m_circles.size());
    }

    /**
     * @brief Number of prevertices K_j on circle j.
     */
    std::size_t prevertexCount(int j) const;

    /**
     * @brief Get the center c_j of circle j.
     */
    Complex getCenter(int j) const;

    /**
     * @brief Get all circle centers, indexed by circle.
     */
    std::vector<Complex> getCenters() const;

    /**
     * @brief Get the radius r_j of circle j.
     */
    double getRadius(int j) const;

    /**
     * @brief Get all circle radii, indexed by circle.
     */
    std::vector<double> getRadii() const;

    /**
     * @brief Get the ordered prevertex angles t(1,j), ..., t(Kj,j) on circle j.
     */
    const std::vector<double>& getPrevertexAngles(int j) const;

    /**
     * @brief Compute the actual prevertex points on circle j: c_j + r_j * exp(i * t(k,j)).
     * @param j Circle index.
     * @return Prevertex points, in prevertex order.
     */
    std::vector<Complex> getPrevertices(int j) const;

    /**
     * @brief Update the center of circle j.
     */
    void setCenter(int j, const Complex& center);

    /**
     * @brief Update the radius of circle j.
     * @throws std::invalid_argument if radius <= 0.
     */
    void setRadius(int j, double radius);

    /**
     * @brief Update the prevertex angles of circle j.
     * @throws std::invalid_argument if the new angle count differs from the current count, or the
     *         angles are not strictly increasing on [t(1,j), t(1,j) + 2*pi).
     */
    void setPrevertexAngles(int j, std::vector<double> angles);

    /**
     * @brief Pack the domain into the flat unconstrained parameter vector Xu used by the MCSC
     *        Newton parameter-problem solve (mirrors circdomain.m's unconstrained()).
     *
     * Packing order (load-bearing -- must match setFromUnconstrained() and any Newton-solve code
     * built on top of this vector):
     *   1. log(r_1), ..., log(r_{m-1})                              -- (m-1) entries, radii of
     *      circles 1..m-1 (circle 0's radius is fixed at 1 by convention and is not packed).
     *   2. Re(c_1), Im(c_1), ..., Re(c_{m-1}), Im(c_{m-1})          -- 2*(m-1) entries, centers of
     *      circles 1..m-1 (circle 0's center is fixed at 0 by convention and is not packed).
     *   3. psi_1,1, ..., psi_{K0-1,0}                                -- (K0-1) entries, unconstrained
     *      angle variables for circle 0 (whose first prevertex angle is fixed at 0 and not packed).
     *   4. For each circle j = 1, ..., m-1, in order:
     *        t(1,j), psi_1,j, ..., psi_{Kj-1,j}                      -- Kj entries: the first
     *        prevertex angle (unconstrained on these circles) followed by its Kj-1 unconstrained
     *        angle variables.
     * Total length: 3*(m-1) - 1 + sum_j K_j = 3*m - 4 + sum_j K_j, matching circdomain.m's Xu.
     */
    Eigen::VectorXd toUnconstrained() const;

    /**
     * @brief Unpack a flat unconstrained parameter vector Xu, updating radii, centers, and
     *        prevertex angles in place (mirrors circdomain.m's set.Xu()). See toUnconstrained()
     *        for the required packing order.
     * @throws std::invalid_argument if Xu's length doesn't match the expected size for the
     *         current circle count and prevertex counts.
     */
    void setFromUnconstrained(const Eigen::VectorXd& Xu);

    bool contains(const Complex& z) const override;

    void transformBoundary(std::function<Complex(const Complex&)> transform) override;

private:
    std::vector<CircleData> m_circles;

    void validateIndex(int j) const;
    static void validateCircleData(const CircleData& circle);

    /**
     * @brief Expected length of the unconstrained parameter vector for the current configuration.
     */
    int expectedUnconstrainedSize() const;
};
