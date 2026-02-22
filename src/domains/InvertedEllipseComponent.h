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

#include "BoundaryComponent.h"

/**
 * @brief Boundary component for the inverted ellipse parameterization (MATLAB binvellip).
 *
 * z(t) = 1 / conj(ellipse(t)) + center, where ellipse(t) is a rotated ellipse
 * with eccentricity controlled by alpha in (0, 1). See design/fornberg/fornmc/binvellip.m.
 */
class InvertedEllipseComponent : public BoundaryComponent
{
public:
    InvertedEllipseComponent(Complex center, double alpha, double rotation = 0.0);

    Complex evaluate(double t) const override;
    Complex evaluateDerivative(double t) const override;
    std::vector<Complex> sample(size_t numPoints) const override;
    double findParameterization(const Complex& z) const override;

    Complex getCenter() const { return m_center; }
    double getAlpha() const { return m_alpha; }
    double getRotation() const { return m_rotation; }

private:
    Complex m_center;
    double m_alpha;
    double m_rotation;
};
