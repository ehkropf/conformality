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

#include "ConformalMapMethod.h"
#include "MCSCReflectionIntegrand.h"
#include "../core/Types.h"
#include "../domains/MCSCCircleDomain.h"
#include "../domains/MCSCPolygonalDomain.h"
#include "../numerics/GaussJacobiQuadrature.h"
#include "../numerics/MCSCContinuationSolver.h"

#include <memory>
#include <optional>
#include <vector>

// Test helper forward declaration
#ifdef TESTING
#include <gtest/gtest_prod.h>
#endif

/**
 * @brief MCSC reflection-method conformal map, unbounded (exterior) case (port of extmap.m's
 *        map_eval, built on the parameter problem solved by #165).
 *
 * compute() solves the unbounded parameter problem (MCSCUnboundedObjectiveFunction, #165) to
 * obtain the circle domain conformally equivalent to the target polygon, then builds the
 * reflection integrand (#163) and Gauss-Jacobi quadrature (#164) used to evaluate the map.
 *
 * map(z) evaluates f(z) by integrating f'(z) along a two-piece path from a known target-polygon
 * vertex image to z (port of extmap.m's integpath: an arc on the source circle from the nearest
 * prevertex to z's angle, then a straight line from the circle boundary out to z), matching
 * mcscmap.m's map_eval + nearestcircle.
 */
class MCSCUnbounded : public ConformalMapMethod
{
#ifdef TESTING
    FRIEND_TEST(MCSCUnboundedTest, NearestCircleFindsContainingCircle);
    FRIEND_TEST(MCSCUnboundedTest, IntegPathSkipsArcAtPrevertex);
    FRIEND_TEST(MCSCUnboundedTest, IntegPathSkipsLineOnBoundary);
#endif

public:
    /**
     * @brief Construct an MCSC unbounded-case method.
     * @param N Reflection truncation level (see MCSCReflectionIntegrand). Defaults to 6, matching
     *        extmapopts.m's reflection-method default.
     * @param ngj Gauss-Jacobi node count per vertex (see GaussJacobiQuadrature).
     * @param continuationOptions Continuation solver tuning for the parameter-problem solve.
     */
    explicit MCSCUnbounded(
        int N = 6,
        int ngj = 12,
        MCSCContinuationSolver::Options continuationOptions = MCSCContinuationSolver::Options()
    );

    ~MCSCUnbounded() override = default;

    /**
     * @brief Solve the parameter problem and build the map-evaluation state.
     * @param map_instance The map to compute; its source domain must be an MCSCCircleDomain
     *        (initial guess) and target domain an MCSCPolygonalDomain.
     * @param target_accuracy Unused -- MCSC's continuation solver has its own internal tolerance
     *        (MCSCContinuationSolver::Options), not a single scalar accuracy target.
     */
    void compute(ConformalMap& map_instance, double target_accuracy = 1e-10) override;

    /**
     * @brief Evaluate f(z): the forward map from the solved circle domain to the target polygon.
     * @param z Point in the source (circle) domain.
     * @return f(z) in the target polygonal domain.
     * @throws std::runtime_error if called before compute().
     */
    Complex map(const Complex& z) const override;

    /**
     * @brief Not implemented -- no MCSC inverse-map issue exists yet.
     * @throws std::runtime_error always.
     */
    Complex inverseMap(const Complex& w) const override;

    /**
     * @brief The solved circle domain, valid after compute().
     * @throws std::runtime_error if called before compute().
     */
    const MCSCCircleDomain& solvedCircleDomain() const;

protected:
    void validateSourceDomain(std::shared_ptr<Domain> domain) const override;
    void validateTargetDomain(std::shared_ptr<Domain> domain) const override;

private:
    int m_N;
    int m_ngj;
    MCSCContinuationSolver::Options m_continuationOptions;

    std::optional<MCSCPolygonalDomain> m_polygon;
    std::optional<MCSCCircleDomain> m_circle;
    std::optional<MCSCReflectionIntegrand> m_integrand;
    std::optional<GaussJacobiQuadrature> m_quadrature;
    Complex m_A{0.0, 0.0};
    std::vector<Complex> m_singularities;   // Flat prevertex list, for the one-half rule
    std::vector<int> m_componentOffset;     // Flat (0-based) vertex index of vertex 0 on component j

    void ensureComputed() const;

    /**
     * @brief Find the circle nearest to z (port of mcscmap.nearestcircle) and the angle of z
     *        relative to that circle's center.
     * @return (circle index, angle in [0, 2*pi)).
     */
    std::pair<int, double> nearestCircle(const Complex& z) const;

    struct ArcSegment
    {
        double left;
        int leftVertex;
        double right;
        int rightVertex;
    };

    struct LineSegment
    {
        Complex left;
        int leftVertex;
    };

    struct IntegrationPath
    {
        std::optional<ArcSegment> arc;
        std::optional<LineSegment> line;
        int nearestVertexFlat;  // Flat index of the nearest prevertex k on circle j
        int circleIndex;        // j
    };

    /**
     * @brief Build the two-piece integration path from the nearest prevertex to z (port of
     *        extmap.m's integpath).
     */
    IntegrationPath integPath(const Complex& z) const;
};
