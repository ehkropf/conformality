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

#include "InvertedEllipseComponent.h"

#include <stdexcept>

InvertedEllipseComponent::InvertedEllipseComponent(Complex center, double alpha, double rotation)
    : m_center(center)
    , m_alpha(alpha)
    , m_rotation(rotation)
{
    if (alpha <= 0.0 || alpha >= 1.0)
    {
        throw std::invalid_argument("InvertedEllipseComponent: alpha must be in (0, 1)");
    }
}

Complex InvertedEllipseComponent::evaluate(double /*t*/) const
{
    throw std::runtime_error("InvertedEllipseComponent::evaluate not implemented");
}

Complex InvertedEllipseComponent::evaluateDerivative(double /*t*/) const
{
    throw std::runtime_error("InvertedEllipseComponent::evaluateDerivative not implemented");
}

std::vector<Complex> InvertedEllipseComponent::sample(size_t /*numPoints*/) const
{
    throw std::runtime_error("InvertedEllipseComponent::sample not implemented");
}

double InvertedEllipseComponent::findParameterization(const Complex& /*z*/) const
{
    throw std::runtime_error("InvertedEllipseComponent::findParameterization not implemented");
}
