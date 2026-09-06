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

#include <complex>
#include <numbers>

// Type aliases for common types used throughout the codebase

using Complex = std::complex<double>;

// Constants used everywhere

constexpr double BOUNDARY_TOLERANCE = 1e-12;

/// Full turn in radians, shared by any code working with angles around a circle.
constexpr double TWO_PI = 2.0 * std::numbers::pi;

// Other

enum class MappingType
{
    INTERIOR_TO_INTERIOR,    // Standard conformal map
    EXTERIOR_TO_INTERIOR,    // External map (most common external case)
    INTERIOR_TO_EXTERIOR,    // Map from bounded to unbounded domain (e.g., unit disk to exterior of unit disk)
    EXTERIOR_TO_EXTERIOR     // Mapping between unbounded domains (e.g., composition of exterior maps)
};

// Forward delcaration to avoid circularity
class Domain;

/// Determines the appropriate mapping type based on the connectivity of source and target domains
MappingType determineMappingType(const Domain& source, const Domain& target);
