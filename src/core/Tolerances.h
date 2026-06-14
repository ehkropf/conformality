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

#include <limits>

/*
 * Centralized numerical tolerances for near-zero checks and division guards.
 *
 * Rationale (see also .claude/references/conformality/tolerances.md):
 *
 * All domains handled by the solver are normalized to the unit disk, so
 * coordinates and the quantities derived from them have O(1) magnitude. With
 * O(1) operands, ABSOLUTE tolerances scaled near machine epsilon (1e-14, 1e-15)
 * are the appropriate choice for geometric-coincidence and pivot guards;
 * converting them to relative tolerances buys nothing and risks regressing the
 * validated thesis-example convergence (thesis2 / thesis4 / thesis5).
 *
 * The one intentionally RELATIVE tolerance is SPLINE_CLOSURE_EPS, because it
 * compares user-supplied input points that are NOT yet normalized, so it is
 * scaled by machine epsilon rather than fixed.
 *
 * Constants are grouped by MEANING, not by value: several share a numeric value
 * today but guard distinct conditions, and should be tuned independently.
 *
 * Do NOT retune any of these without re-running the thesis regression tests --
 * a change in their numerical output indicates the tolerance mattered.
 *
 * Related: BOUNDARY_TOLERANCE (1e-12) in Types.h is the domain-containment
 * tolerance and lives there with the other domain-level definitions.
 */

/// Two points/centers are treated as coincident, or a complex magnitude treated
/// as zero, on the unit-disk-normalized domain (Laurent-series denominators,
/// hole-center coincidence, Sherman-Morrison singularity).
constexpr double GEOMETRIC_COINCIDENCE_EPS = 1e-14;

/// Guards the per-point tangent-magnitude divisor in the Fornberg radius
/// (Newton) update. Distinct in meaning from GEOMETRIC_COINCIDENCE_EPS even
/// though equal today: it guards a different physical quantity (the update
/// scaling), and the two should be tuned independently.
constexpr double NEWTON_UPDATE_SCALE_EPS = 1e-14;

/// Pivot / segment-length / derivative guards that protect a division
/// (ray-cast segment height, spline chord and tridiagonal pivots, Newton
/// derivative breakdown).
constexpr double PIVOT_EPS = 1e-15;

/// Endpoint-coincidence (closed-curve) test on raw, un-normalized user input.
/// Intentionally RELATIVE: scaled by machine epsilon rather than a fixed value.
constexpr double SPLINE_CLOSURE_EPS = 100.0 * std::numeric_limits<double>::epsilon();

/// Grid geometry on the unit-disk-scaled domain: degenerate-radius clamp,
/// full-circle detection, and point deduplication.
constexpr double GRID_GEOMETRY_EPS = 1e-10;

/// Step size for central/forward finite-difference derivatives and normals
/// (a deliberate ~sqrt(machine epsilon) choice).
constexpr double FINITE_DIFFERENCE_STEP = 1e-6;
