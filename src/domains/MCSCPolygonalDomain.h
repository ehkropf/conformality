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
#include "Domain.h"

#include <functional>
#include <vector>

/**
 * @brief Target polygonal domain for the MCSC method (reflection method, Phase 5)
 *
 * This is the polygonal analog of the canonical (source) domain for MCSC: a thin,
 * effectively immutable data holder for the m boundary-component vertex lists and
 * their interior angles, mirroring the MATLAB reference's `polygon.m`/`polygons.m`.
 *
 * Unlike `SimplyConnectedDomain`/`MultiplyConnectedDomain`, this type does not use
 * `Boundary`/`BoundaryComponent` at all -- MCSC's reflection method enforces
 * conditions at discrete prevertex/collocation points and has no continuous
 * parametrization or FFT/dense-sampling step. No unknowns/mutable Newton-solve
 * state live here; the Newton solve's unknowns are entirely circle-side
 * (`MCSCCircleDomain`). This class is pure data plus validation.
 *
 * Orientation convention: each component's vertex list is independently checked and
 * (if needed) re-oriented to be counterclockwise as a standalone simple closed curve
 * (mirrors `polygon.m`'s per-polygon `calc_angles`, applied identically to the outer
 * boundary and every hole). This class does not apply the additional outer/inner sign
 * asymmetry the MATLAB reference layers on top in `intpolys.m`/`extpolys.m` (negating
 * `beta` for the outer component only) -- that is a mapping-formulation concern for a
 * later issue, not a property of the vertex/angle data stored here.
 */
class MCSCPolygonalDomain : public Domain
{
private:
    std::vector<std::vector<Complex>> m_vertices;  // Per-component vertex lists (each CCW; see above)
    std::vector<std::vector<double>> m_alpha;      // Per-component interior angles (as multiples of pi)

public:
    /**
     * @brief Construct a target polygonal domain from per-component vertex lists
     * @param vertices Per-boundary-component ordered vertex lists, each independently a
     *        simple closed curve (orientation is normalized to counterclockwise; see
     *        class-level orientation convention notes)
     * @param isUnboundedDomain Whether this is an unbounded (exterior) target domain
     *
     * Interior angles are computed and validated for each component (mirrors
     * `polygon.m`'s `calc_angles`): each component's vertex list is checked against
     * `sum(1-alpha) == 2` (up to a `100*sqrt(n)*eps` tolerance) and re-oriented
     * (vertices and angles reversed) if the raw angle sum is negative, matching
     * `polygon.m`'s orientation-fixup convention.
     */
    explicit MCSCPolygonalDomain(
        const std::vector<std::vector<Complex>>& vertices,
        bool isUnboundedDomain = false
    );

    ~MCSCPolygonalDomain() override = default;

    /**
     * @brief Get the vertex lists for all boundary components
     * @return Vector of per-component vertex lists
     */
    const std::vector<std::vector<Complex>>& getVertices() const
    {
        return m_vertices;
    }

    /**
     * @brief Get the vertex list for a single boundary component
     * @param componentIndex Index of the boundary component
     * @return Vertex list for the component
     */
    const std::vector<Complex>& getVertices(int componentIndex) const;

    /**
     * @brief Get the interior angles (as multiples of pi) for all boundary components
     * @return Vector of per-component interior-angle lists
     */
    const std::vector<std::vector<double>>& getAlpha() const
    {
        return m_alpha;
    }

    /**
     * @brief Get the interior angles (as multiples of pi) for a single boundary component
     * @param componentIndex Index of the boundary component
     * @return Interior-angle list for the component
     */
    const std::vector<double>& getAlpha(int componentIndex) const;

    /**
     * @brief Get the vertex count for a single boundary component
     * @param componentIndex Index of the boundary component
     * @return Number of vertices in the component
     */
    int vertexCount(int componentIndex) const;

    bool contains(const Complex& z) const override;

    /**
     * @brief Transform all vertices using a function, recomputing angles afterward
     * @param transform Function mapping complex points to complex points
     * @throws std::invalid_argument if any transformed component fails polygon validation
     *
     * Applies transform to every vertex in every component, then recomputes and
     * validates interior angles from the transformed vertex lists, re-orienting any
     * component whose transformed vertex list is wound clockwise (e.g. an
     * orientation-reversing transform like conjugation) -- the same per-component
     * normalization applied at construction time. This method is atomic: if any
     * component fails validation, the domain is left completely unmodified (as if
     * the call had not been made) rather than partially transformed.
     */
    void transformBoundary(std::function<Complex(const Complex&)> transform) override;

private:
    /**
     * @brief Compute and validate interior angles for a single component's vertex list
     * @param vertices Ordered vertex list for one boundary component
     * @return Interior angles (as multiples of pi), one per vertex
     * @throws std::invalid_argument if the vertex list does not form a valid polygon
     *
     * Port of `polygon.m`'s static `calc_angles`. Computes the exterior turning angle
     * at each vertex, checks that the total angle sum matches a valid simple polygon
     * (`sum(1-alpha) == 2`, up to tolerance), and reverses (vertices + angles) if the
     * raw sum indicates the wrong orientation.
     */
    static std::vector<double> calcAngles(std::vector<Complex>& vertices);
};
