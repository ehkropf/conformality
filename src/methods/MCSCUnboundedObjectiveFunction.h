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

#include "../domains/MCSCCircleDomain.h"
#include "../domains/MCSCPolygonalDomain.h"
#include "../numerics/GaussJacobiQuadrature.h"
#include "../numerics/MCSCContinuationSolver.h"
#include "MCSCReflectionIntegrand.h"

#include <Eigen/Dense>

/**
 * @brief MCSC unbounded-case parameter-problem residual function (dissertation Sec 3.1.1; port of
 *        extobjfun.m's fun_eval).
 *
 * Given a target (unbounded) polygonal domain and the unconstrained-parameter vector Xu for a
 * candidate circle domain, evaluate() computes the side-length, position, and orientation
 * residuals that vanish exactly when the circle domain is the true conformally-equivalent
 * pre-image. This is the F(Xu) callable driving MCSCContinuationSolver's homotopy solve -- the
 * direct analog of FornbergMC::formSystem/solveSystem/newtonUpdate for the MCSC method.
 *
 * The residual ordering below must match MCSCCircleDomain::toUnconstrained()'s packing convention
 * exactly -- the two are not independently reorderable, since the continuation solver treats
 * MCSCCircleDomain's Xu as an opaque R^n vector with no knowledge of what each entry means:
 *   1. For circles j = 1, ..., m-1: orientation + first-side-length condition (2 real equations
 *      each -- Re/Im of the complex condition SL(0,j) - (w(1,j) - w(0,j)) = 0). Circle 0's first
 *      side is not constrained here -- it is fixed by definition of the scaling constant A.
 *   2. For circles j = 1, ..., m-1: position condition relative to circle 0 (2 real equations
 *      each -- Re/Im of T(j-1) - (w(0,j) - w(0,0)) = 0).
 *   3. Remaining side-length conditions (magnitude only, 1 real equation each), for every circle
 *      including circle 0, skipping each circle's first side (already covered by step 1 for
 *      j >= 1, or fixed by A for j = 0).
 *
 * Deviation from extobjfun.m: the MATLAB reference's position-condition loop reuses a stale loop
 * variable ('j', left over at value m from the preceding side-length loop) when computing the
 * flat prevertex index passed to lineq's rpvn argument, so every position-condition integral is
 * (except when m == 2) evaluated with the wrong right-endpoint vertex index. This port computes
 * the correct per-j flat index instead -- see MCSCUnboundedObjectiveFunction.cpp for detail.
 */
class MCSCUnboundedObjectiveFunction
{
public:
    /**
     * @brief Construct the objective function around a working circle domain.
     * @param polygon Target (unbounded) polygonal domain -- fixed for the life of the solve; must
     *        outlive this object.
     * @param initialCircle Initial circle-domain configuration -- copied into an internal working
     *        domain that evaluate()/solve() mutate. Connectivity and per-circle prevertex counts
     *        must match polygon's connectivity and per-component vertex counts.
     * @param N Reflection truncation level passed to MCSCReflectionIntegrand. Defaults to 6,
     *        matching extmapopts.m's fpclasstable default for the reflection method ('refl' ->
     *        fpextrefl, N=6) -- not to be confused with the unrelated fast/least-squares method's
     *        default of 25 ('fpls' -> fpextfpls), which is out of Phase 5's scope. Since the
     *        reflection sequence grows as sum_{level=0}^{N} (m-1)^level per circle, a much larger
     *        N becomes intractable quickly for m > 2 (e.g. N=25, m=3 is ~2*10^8 reflected circles
     *        per original circle).
     * @param ngj Gauss-Jacobi node count per vertex (see GaussJacobiQuadrature).
     * @throws std::invalid_argument under the same conditions as MCSCReflectionIntegrand's
     *         constructor (connectivity/prevertex-count mismatch between polygon and
     *         initialCircle).
     */
    MCSCUnboundedObjectiveFunction(
        const MCSCPolygonalDomain& polygon, const MCSCCircleDomain& initialCircle, int N = 6, int ngj = 12
    );

    /**
     * @brief Evaluate F(Xu): the side-length/position/orientation residual vector.
     *
     * Mutates the working circle domain (see workingCircleDomain()) to the configuration encoded
     * by Xu before evaluating -- matching extobjfun.m's `O.C.Xu = Xu` side effect.
     *
     * @param Xu Unconstrained parameter vector (see MCSCCircleDomain::toUnconstrained()).
     * @return Residual vector, same length as Xu.
     * @throws std::invalid_argument if Xu's length doesn't match the working circle domain's
     *         expected unconstrained-vector length (propagated from
     *         MCSCCircleDomain::setFromUnconstrained).
     */
    Eigen::VectorXd evaluate(const Eigen::VectorXd& Xu);

    /**
     * @brief Solve the parameter problem via numerical continuation (dissertation Sec 3.3),
     *        starting from the working circle domain's current configuration.
     * @param options Continuation solver tuning (defaults match MCSCContinuationSolver::Options).
     * @return The solved circle domain.
     * @throws MCSCContinuationSolver::ConvergenceError if continuation fails to converge.
     */
    MCSCCircleDomain solve(MCSCContinuationSolver::Options options = MCSCContinuationSolver::Options());

    /**
     * @brief The working circle domain, as last set by evaluate(). Provided so callers (e.g.
     *        tests, or solve()'s own caller) can inspect the domain corresponding to a given Xu
     *        without re-solving setFromUnconstrained() themselves.
     */
    const MCSCCircleDomain& workingCircleDomain() const
    {
        return m_circle;
    }

private:
    const MCSCPolygonalDomain& m_polygon;
    int m_N;
    MCSCCircleDomain m_circle;
    MCSCReflectionIntegrand m_integrand;
    GaussJacobiQuadrature m_quadrature;

    static GaussJacobiQuadrature makeQuadrature(const MCSCPolygonalDomain& polygon, int ngj);
};
