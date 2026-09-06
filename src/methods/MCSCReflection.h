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

#include <vector>

/**
 * @brief Circle-reflection geometry shared by the MCSC reflection method, bounded and unbounded
 *        alike (dissertation Sec 2.1; port of fprefl.m's static methods and reflectzsmi()).
 *
 * These are pure functions on circle data -- they know nothing about polygons, prevertex
 * correspondence, or which of the bounded/unbounded cases is in play. Both fpextrefl.m and
 * fpintrefl.m in the MATLAB reference call this exact same reflection machinery; the bounded case
 * (#167) additionally reflects its interior circles to the exterior before calling
 * reflectCircleSequence(), a step that belongs to that case, not to this file.
 */
namespace mcsc
{

/**
 * @brief Reflect circle (ci, ri) through circle (c, r) to get circle (co, ro).
 *
 * Port of fprefl.m's centrad(). Follows directly from the classical circle-inversion formula: a
 * circle centered at ci with radius ri, reflected through the circle (c, r), maps to another
 * circle (not, in general, centered at the image of ci).
 *
 * @param c Center of the reflecting circle.
 * @param r Radius of the reflecting circle (must be positive).
 * @param ci Center of the circle being reflected.
 * @param ri Radius of the circle being reflected (must be positive).
 * @return Pair (co, ro): center and radius of the reflected circle.
 * @throws std::invalid_argument if r or ri is non-positive, or if (ci, ri) is centered exactly at
 *         c with |ci - c| == ri (division by zero in the reflection formula).
 */
std::pair<Complex, double> reflectCircle(const Complex& c, double r, const Complex& ci, double ri);

/**
 * @brief Reflect point z through circle (c, r) to get the reflected point.
 *
 * Port of fprefl.m's reflect(): zo = c + r^2 / conj(z - c).
 *
 * @param c Center of the reflecting circle.
 * @param r Radius of the reflecting circle (must be positive).
 * @param z Point to reflect.
 * @return Reflected point.
 * @throws std::invalid_argument if r is non-positive, or if z == c (division by zero).
 */
Complex reflectPoint(const Complex& c, double r, const Complex& z);

/**
 * @brief One reflected copy of the original m-circle configuration, at some reflection level.
 *
 * Mirrors one "column" of the MATLAB arrays cnu(nu,j), rnu(nu,j), znu(k,nu,j), snu(nu,j) for a
 * fixed circle index j and reflection index nu (nu = 0 is the original, unreflected circle).
 */
struct ReflectedCircle
{
    Complex center;                    ///< c_nu
    double radius{0.0};                ///< r_nu
    std::vector<Complex> prevertices;  ///< z_{k,nu} for this circle's own prevertex count
    Complex outerImage;                ///< s_nu -- the reflected image of the outer circle's center
    int lastReflectedThrough{-1};      ///< jlr(nu,j): index of the circle this copy was last reflected through
};

/**
 * @brief Build the truncated reflection sequence for every circle out to level N (dissertation
 *        Sec 2.1; port of fprefl.m's reflectzsmi(), method-of-images variant).
 *
 * For each of the m original circles j, and its own prevertices, repeatedly reflects through
 * every other circle in turn -- never immediately re-reflecting through the circle a copy was
 * last produced from -- out to N levels. Level 0 is the original, unreflected circle.
 *
 * @param centers Original circle centers c_1, ..., c_m (size m).
 * @param radii Original circle radii r_1, ..., r_m (size m, all positive).
 * @param prevertices Original per-circle prevertex points z_{k,j} (size m; prevertices[j].size()
 *        may vary per circle).
 * @param outerCenters Original per-circle "outer image" points to also reflect alongside the
 *        prevertices (size m) -- for the unbounded case this is the circle's own center (see
 *        fpextrefl.build_reflections's cn passed as both cn and sn); the bounded case (#167)
 *        supplies a different value here after its own pre-reflection step.
 * @param N Truncation level (number of reflection levels beyond the original circles, N >= 0).
 * @return Per-circle list of ReflectedCircle, one entry per reflection level (size
 *         1 + m + m*(m-1) + ... up to level N, i.e. sum_{level=0}^{N} (m-1)^level for level >= 1
 *         and 1 for level 0, matching reflectzsmi()'s (m-1)^level growth).
 * @throws std::invalid_argument if centers/radii/prevertices/outerCenters have inconsistent
 *         sizes, if any radius is non-positive, if there are fewer than 2 circles (a reflection
 *         method needs at least one other circle to reflect through), or if N is negative.
 */
std::vector<std::vector<ReflectedCircle>> reflectCircleSequence(
    const std::vector<Complex>& centers,
    const std::vector<double>& radii,
    const std::vector<std::vector<Complex>>& prevertices,
    const std::vector<Complex>& outerCenters,
    int N
);

} // namespace mcsc
