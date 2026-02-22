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

#include "BoundaryHelpers.h"

#include <cmath>
#include <stdexcept>

namespace conformality::examples
{

std::shared_ptr<Boundary> createCircularBoundary(Complex center, double radius)
{
    if (radius <= 0.0)
    {
        throw std::invalid_argument("createCircularBoundary: radius must be positive");
    }
    auto component = std::make_shared<AnalyticBoundaryComponent>(
        [center, radius](double theta) {
            return center + radius * Complex(std::cos(theta), std::sin(theta));
        },
        [radius](double theta) {
            return radius * Complex(-std::sin(theta), std::cos(theta));
        }
    );
    return std::make_shared<Boundary>(component);
}

std::shared_ptr<Boundary> createEllipseBoundary(
    Complex center, double semi_major, double semi_minor, double rotation)
{
    if (semi_major <= 0.0 || semi_minor <= 0.0)
    {
        throw std::invalid_argument("createEllipseBoundary: semi-axes must be positive");
    }
    auto component = std::make_shared<AnalyticBoundaryComponent>(
        [center, semi_major, semi_minor, rotation](double t) {
            double cosR = std::cos(rotation);
            double sinR = std::sin(rotation);
            double x = semi_major * std::cos(t);
            double y = semi_minor * std::sin(t);
            return center + Complex(x * cosR - y * sinR, x * sinR + y * cosR);
        },
        [semi_major, semi_minor, rotation](double t) {
            double cosR = std::cos(rotation);
            double sinR = std::sin(rotation);
            double dx = -semi_major * std::sin(t);
            double dy = semi_minor * std::cos(t);
            return Complex(dx * cosR - dy * sinR, dx * sinR + dy * cosR);
        }
    );
    return std::make_shared<Boundary>(component);
}

std::shared_ptr<Boundary> createInvertedEllipseBoundary(
    Complex center, double alpha, double rotation)
{
    auto component = std::make_shared<InvertedEllipseComponent>(center, alpha, rotation);
    return std::make_shared<Boundary>(component);
}

} // namespace conformality::examples
