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
#include "../numerics/RootFinder.h"

#include <cmath>
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

Complex InvertedEllipseComponent::evaluate(double t) const
{
    double cosS = std::cos(t);
    double sinS = std::sin(t);
    double cosR = std::cos(m_rotation);
    double sinR = std::sin(m_rotation);
    double a_cosR = m_alpha * cosR;
    double a_sinR = m_alpha * sinR;

    Complex conj_ell(a_cosR * cosS - sinR * sinS,
                     -(cosR * sinS + a_sinR * cosS));

    return 1.0 / conj_ell + m_center;
}

Complex InvertedEllipseComponent::evaluateDerivative(double t) const
{
    double cosS = std::cos(t);
    double sinS = std::sin(t);
    double cosR = std::cos(m_rotation);
    double sinR = std::sin(m_rotation);
    double a_cosR = m_alpha * cosR;
    double a_sinR = m_alpha * sinR;

    Complex conj_ell(a_cosR * cosS - sinR * sinS,
                     -(cosR * sinS + a_sinR * cosS));

    Complex numerator(a_cosR * sinS + sinR * cosS,
                      cosR * cosS - a_sinR * sinS);

    return numerator / (conj_ell * conj_ell);
}

std::vector<Complex> InvertedEllipseComponent::sample(size_t numPoints) const
{
    std::vector<Complex> samples;
    samples.reserve(numPoints);
    for (size_t i = 0; i < numPoints; ++i)
    {
        double t = 2.0 * M_PI * i / numPoints;
        samples.push_back(evaluate(t));
    }
    return samples;
}

double InvertedEllipseComponent::findParameterization(const Complex& z) const
{
    auto objective = [this, &z](double t) -> double
    {
        Complex diff = this->evaluate(t) - z;
        return std::norm(diff);
    };

    try
    {
        double result = RootFinder::ternarySearch(objective, 0.0, 2.0 * M_PI, 1e-12);
        while (result < 0.0) result += 2.0 * M_PI;
        while (result >= 2.0 * M_PI) result -= 2.0 * M_PI;
        return result;
    }
    catch (const RootFinder::ConvergenceError&)
    {
        if (p_statusManager)
        {
            p_statusManager->reportWarning("InvertedEllipseComponent",
                                           "Root finding failed to converge in findParameterization");
        }
        throw;
    }
}
