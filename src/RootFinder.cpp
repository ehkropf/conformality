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

#include "RootFinder.h"

#include <cmath>

double RootFinder::ternarySearch(std::function<double(double)> objective,
                                 double low, double high,
                                 double tolerance,
                                 int maxIterations)
{
    if (low >= high)
    {
        throw std::invalid_argument("Invalid interval: low must be less than high");
    }

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        double range = std::abs(high - low);
        double scale = std::max(std::abs(high), std::abs(low));
        if (range < tolerance * std::max(scale, 1.0))
        {
            return (low + high) / 2.0;
        }

        double delta = (high - low) / 3.0;
        double mid1 = low + delta;
        double mid2 = high - delta;

        double f1 = objective(mid1);
        double f2 = objective(mid2);

        if (f1 > f2)
        {
            low = mid1;
        }
        else
        {
            high = mid2;
        }
    }

    throw ConvergenceError("Ternary search failed to converge within " +
                           std::to_string(maxIterations) + " iterations");
}

