/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
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

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "BoundaryComponent.hpp"
#include <memory>
#include <limits>
#include <stdexcept>

/**
 * @brief Abstract base class for boundaries
 *
 * A boundary consists of one or more boundary components
 */
class Boundary
{
protected:
    std::vector<std::shared_ptr<BoundaryComponent>> components;

public:
    virtual ~Boundary() = default;

    /**
     * @brief Evaluate the boundary at parameter t
     * @param t Parameter value
     * @return ComplexDouble point on the boundary
     */
    virtual ComplexDouble evaluate(double t) const = 0;

    /**
     * @brief Evaluate the derivative at parameter t
     * @param t Parameter value
     * @return ComplexDouble derivative
     */
    virtual ComplexDouble evaluateDerivative(double t) const = 0;

    /**
     * @brief Sample points along the boundary
     * @param numPoints Number of points to sample
     * @return Vector of complex points on the boundary
     */
    virtual std::vector<ComplexDouble> sample(int numPoints) const = 0;

    /**
     * @brief Find the parameter value for a point on the boundary
     * @param z ComplexDouble point on or near the boundary
     * @return Parameter value t such that evaluate(t) ≈ z
     */
    virtual double findParameterization(const ComplexDouble& z) const = 0;

    /**
     * @brief Add a boundary component
     * @param component Boundary component to add
     */
    void addComponent(std::shared_ptr<BoundaryComponent> component)
    {
        component->setIndex(components.size());
        components.push_back(component);
    }

    /**
     * @brief Get the number of boundary components
     * @return Number of components
     */
    size_t getNumComponents() const
    {
        return components.size();
    }

    /**
     * @brief Get a specific boundary component
     * @param index Component index
     * @return Boundary component
     */
    const BoundaryComponent& getComponent(size_t index) const
    {
        if (index >= components.size())
        {
            throw std::out_of_range("Boundary component index out of range");
        }
        return *components[index];
    }

    /**
     * @brief Get all boundary components
     * @return Vector of boundary components
     */
    const std::vector<std::shared_ptr<BoundaryComponent>>& getComponents() const
    {
        return components;
    }
};

/**
 * @brief Simple boundary consisting of a single component
 */
class SimpleBoundary : public Boundary
{
public:
    /**
     * @brief Construct a new Simple Boundary
     * @param component Single boundary component
     */
    SimpleBoundary(std::shared_ptr<BoundaryComponent> component)
    {
        addComponent(component);
    }

    ComplexDouble evaluate(double t) const override
    {
        return components[0]->evaluate(t);
    }

    ComplexDouble evaluateDerivative(double t) const override
    {
        return components[0]->evaluateDerivative(t);
    }

    std::vector<ComplexDouble> sample(int numPoints) const override
    {
        return components[0]->sample(numPoints);
    }

    double findParameterization(const ComplexDouble& z) const override
    {
        return components[0]->findParameterization(z);
    }
};

/**
 * @brief Composite boundary consisting of multiple components
 */
class CompositeBoundary : public Boundary
{
private:
    std::vector<double> parameterRanges;  // End of parameter range for each component

public:
    CompositeBoundary() = default;

    /**
     * @brief Construct a new Composite Boundary
     * @param components Vector of boundary components
     */
    CompositeBoundary(std::vector<std::shared_ptr<BoundaryComponent>> boundaryComponents)
    {
        for (auto& component : boundaryComponents)
        {
            addComponent(component);
        }
        updateParameterRanges();
    }

    void addComponent(std::shared_ptr<BoundaryComponent> component)
    {
        Boundary::addComponent(component);
        updateParameterRanges();
    }

    ComplexDouble evaluate(double t) const override
    {
        double normalizedT = std::fmod(t, 2.0 * M_PI);
        if (normalizedT < 0) normalizedT += 2.0 * M_PI;

        // Find which component the parameter corresponds to
        size_t componentIdx = 0;
        double componentT = normalizedT;

        for (size_t i = 0; i < parameterRanges.size(); ++i)
        {
            if (normalizedT < parameterRanges[i])
            {
                componentIdx = i;
                break;
            }
            componentT -= parameterRanges[i];
        }

        // Scale parameter to component's range [0, 2π)
        double scaledT = componentT * 2.0 * M_PI / parameterRanges[componentIdx];
        return components[componentIdx]->evaluate(scaledT);
    }

    ComplexDouble evaluateDerivative(double t) const override
    {
        double normalizedT = std::fmod(t, 2.0 * M_PI);
        if (normalizedT < 0) normalizedT += 2.0 * M_PI;

        // Find which component the parameter corresponds to
        size_t componentIdx = 0;
        double componentT = normalizedT;

        for (size_t i = 0; i < parameterRanges.size(); ++i)
        {
            if (normalizedT < parameterRanges[i])
            {
                componentIdx = i;
                break;
            }
            componentT -= parameterRanges[i];
        }

        // Scale parameter to component's range [0, 2π)
        double scaledT = componentT * 2.0 * M_PI / parameterRanges[componentIdx];
        return components[componentIdx]->evaluateDerivative(scaledT);
    }

    std::vector<ComplexDouble> sample(int numPoints) const override
    {
        if (components.empty() || numPoints <= 0)
        {
            return {};
        }

        std::vector<ComplexDouble> samples;
        samples.reserve(numPoints);

        // Distribute points proportionally to parameter ranges
        for (size_t i = 0; i < components.size(); ++i)
        {
            int componentPoints = static_cast<int>(numPoints * parameterRanges[i] / (2.0 * M_PI));
            if (componentPoints <= 0)
            {
                componentPoints = 1;  // Ensure at least one point per component
            }

            auto componentSamples = components[i]->sample(componentPoints);
            samples.insert(samples.end(), componentSamples.begin(), componentSamples.end());
        }

        // Adjust to exactly numPoints if needed
        if (samples.size() > static_cast<size_t>(numPoints))
        {
            samples.resize(numPoints);
        }
        else if (samples.size() < static_cast<size_t>(numPoints))
        {
            // Add more points to the first component if needed
            int remaining = numPoints - samples.size();
            auto extraSamples = components[0]->sample(remaining);
            samples.insert(samples.end(), extraSamples.begin(), extraSamples.end());
        }

        return samples;
    }

    double findParameterization(const ComplexDouble& z) const override
    {
        if (components.empty())
        {
            return 0.0;
        }

        // Find closest component and its parameter
        double minDist = std::numeric_limits<double>::max();
        double bestParam = 0.0;

        for (size_t i = 0; i < components.size(); ++i)
        {
            double t = components[i]->findParameterization(z);
            ComplexDouble point = components[i]->evaluate(t);
            double dist = std::norm(z.getValue() - point.getValue());

            if (dist < minDist)
            {
                minDist = dist;
                // Map local parameter to global parameter
                double globalParam = 0.0;
                for (size_t j = 0; j < i; ++j)
                {
                    globalParam += parameterRanges[j];
                }
                globalParam += t * parameterRanges[i] / (2.0 * M_PI);
                bestParam = globalParam;
            }
        }

        return bestParam;
    }

private:
    void updateParameterRanges()
    {
        parameterRanges.resize(components.size());

        // Simple approach: divide parameter space equally
        double paramPerComponent = 2.0 * M_PI / components.size();

        for (size_t i = 0; i < components.size(); ++i)
        {
            parameterRanges[i] = paramPerComponent;
        }

        // Convert to cumulative ranges
        for (size_t i = 1; i < parameterRanges.size(); ++i)
        {
            parameterRanges[i] += parameterRanges[i-1];
        }
    }
};

#endif // BOUNDARY_HPP
