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

#include <Eigen/Dense>
#include <cassert>
#include <complex>

/**
 * @brief Encapsulates conformal moduli for multiply-connected domains
 *
 * The conformal moduli are the centers (c) and radii (rho) of the circular
 * holes in the canonical domain representation. For an m-connected domain:
 * - c has m-1 entries (centers of inner circles)
 * - rho has m-1 entries (radii of inner circles)
 */
struct ConformalModuli
{
    Eigen::VectorXcd c;    // Centers of inner circles (m-1 values)
    Eigen::VectorXd rho;   // Radii of inner circles (m-1 values)

    /// Access center of inner circle L (0-based index)
    std::complex<double> center(int L) const
    {
        assert(L >= 0 && L < c.size());
        return c(L);
    }

    /// Access radius of inner circle L (0-based index)
    double radius(int L) const
    {
        assert(L >= 0 && L < rho.size());
        return rho(L);
    }

    /// Number of inner circles
    int numInnerCircles() const { return static_cast<int>(rho.size()); }
};
