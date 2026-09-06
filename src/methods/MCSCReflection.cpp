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

#include "MCSCReflection.h"

#include <stdexcept>

namespace mcsc
{

std::pair<Complex, double> reflectCircle(const Complex& c, double r, const Complex& ci, double ri)
{
    if (r <= 0.0 || ri <= 0.0)
    {
        throw std::invalid_argument("mcsc::reflectCircle: circle radii must be positive");
    }

    const Complex diff = ci - c;
    const double denom = std::norm(diff) - ri * ri;  // |ci - c|^2 - ri^2
    if (denom == 0.0)
    {
        throw std::invalid_argument("mcsc::reflectCircle: circle to reflect passes through the reflecting center");
    }

    const Complex co = c + r * r * diff / denom;
    const double ro = r * r * ri / std::abs(denom);
    return {co, ro};
}

Complex reflectPoint(const Complex& c, double r, const Complex& z)
{
    if (r <= 0.0)
    {
        throw std::invalid_argument("mcsc::reflectPoint: reflecting circle radius must be positive");
    }
    if (z == c)
    {
        throw std::invalid_argument("mcsc::reflectPoint: cannot reflect the reflecting circle's own center");
    }

    return c + r * r / std::conj(z - c);
}

namespace
{

void validateReflectCircleSequenceInputs(
    const std::vector<Complex>& centers,
    const std::vector<double>& radii,
    const std::vector<std::vector<Complex>>& prevertices,
    const std::vector<Complex>& outerCenters,
    int N
)
{
    const std::size_t m = centers.size();
    if (m < 2)
    {
        throw std::invalid_argument("mcsc::reflectCircleSequence: at least 2 circles are required");
    }
    if (radii.size() != m || prevertices.size() != m || outerCenters.size() != m)
    {
        throw std::invalid_argument("mcsc::reflectCircleSequence: centers/radii/prevertices/outerCenters size mismatch");
    }
    for (double r : radii)
    {
        if (r <= 0.0)
        {
            throw std::invalid_argument("mcsc::reflectCircleSequence: circle radii must be positive");
        }
    }
    if (N < 0)
    {
        throw std::invalid_argument("mcsc::reflectCircleSequence: truncation level N must be non-negative");
    }
}

} // namespace

std::vector<std::vector<ReflectedCircle>> reflectCircleSequence(
    const std::vector<Complex>& centers,
    const std::vector<double>& radii,
    const std::vector<std::vector<Complex>>& prevertices,
    const std::vector<Complex>& outerCenters,
    int N
)
{
    validateReflectCircleSequenceInputs(centers, radii, prevertices, outerCenters, N);

    const int m = static_cast<int>(centers.size());
    std::vector<std::vector<ReflectedCircle>> result(m);

    for (int j = 0; j < m; ++j)
    {
        ReflectedCircle level0;
        level0.center = centers[j];
        level0.radius = radii[j];
        level0.prevertices = prevertices[j];
        level0.outerImage = outerCenters[j];
        level0.lastReflectedThrough = j;
        result[j].push_back(level0);

        // Breadth-first expansion: at each level, every copy from the previous level spawns
        // m-1 children, one per circle j1 != the circle it was last reflected through (a
        // reflection sequence never immediately re-reflects through the same circle, matching
        // reflectzsmi()'s jlr bookkeeping).
        std::size_t levelStart = 0;
        std::size_t levelEnd = result[j].size();
        for (int level = 1; level <= N; ++level)
        {
            std::vector<ReflectedCircle> newChildren;
            for (std::size_t nu = levelStart; nu < levelEnd; ++nu)
            {
                // Copy (not reference) the parent: result[j] is appended to below (via
                // newChildren, merged in after this loop), and taking a reference into a vector
                // that gets reallocated during the same loop would dangle.
                const ReflectedCircle parent = result[j][nu];
                for (int j1 = 0; j1 < m; ++j1)
                {
                    if (j1 == parent.lastReflectedThrough)
                    {
                        continue;
                    }

                    ReflectedCircle child;
                    const auto [co, ro] = reflectCircle(centers[j1], radii[j1], parent.center, parent.radius);
                    child.center = co;
                    child.radius = ro;
                    child.outerImage = reflectPoint(centers[j1], radii[j1], parent.outerImage);
                    child.lastReflectedThrough = j1;

                    child.prevertices.reserve(parent.prevertices.size());
                    for (const Complex& z : parent.prevertices)
                    {
                        child.prevertices.push_back(reflectPoint(centers[j1], radii[j1], z));
                    }

                    newChildren.push_back(std::move(child));
                }
            }
            levelStart = levelEnd;
            result[j].insert(result[j].end(), newChildren.begin(), newChildren.end());
            levelEnd = result[j].size();
        }
    }

    return result;
}

} // namespace mcsc
