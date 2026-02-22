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

#include <stdexcept>

namespace conformality::examples
{

std::shared_ptr<Boundary> createCircularBoundary(Complex /*center*/, double /*radius*/)
{
    throw std::runtime_error("createCircularBoundary not implemented");
}

std::shared_ptr<Boundary> createEllipseBoundary(
    Complex /*center*/, double /*semi_major*/, double /*semi_minor*/, double /*rotation*/)
{
    throw std::runtime_error("createEllipseBoundary not implemented");
}

std::shared_ptr<Boundary> createInvertedEllipseBoundary(
    Complex /*center*/, double /*alpha*/, double /*rotation*/)
{
    throw std::runtime_error("createInvertedEllipseBoundary not implemented");
}

} // namespace conformality::examples
