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

#ifndef GRID_HPP
#define GRID_HPP

#include "../core/Types.h"
#include <vector>
#include <functional>
#include <unordered_map>
#include <string>

/**
 * @brief Enumeration of available grid types
 */
enum class GridType
{
    POLAR,
    CARTESIAN,
    PARAMETRIC
};

/**
 * @brief Class for representing grid points in a complex plane
 *
 * This class can handle various types of grids (polar, Cartesian, etc.) and
 * supports transformations between different grids via conformal mappings.
 */
class Grid
{
private:
    GridType type;
    std::vector<Complex> points;
    std::vector<std::vector<size_t>> lines;

    // Parameters for grid generation
    double xMin, xMax, yMin, yMax;     // Cartesian grid bounds
    double rMin, rMax;                 // Polar grid radial bounds
    double angleMin, angleMax;         // Polar grid angular bounds

    int numRadialLines;                // For polar grids
    int numAngularLines;               // For polar grids
    int numHorizontalLines;            // For Cartesian grids
    int numVerticalLines;              // For Cartesian grids

    bool includeOrigin;                // Whether to include origin in polar grids

    // Parametric grid properties
    std::function<Complex(double, double)> parameterization;
    double uMin, uMax, vMin, vMax;     // Parametric grid bounds
    int numULines, numVLines;          // Parametric grid line counts

    /**
     * @brief Validate a grid line count, rejecting degenerate single-line grids
     *
     * Grid spacing is computed as (max - min) / (count - 1), so a count of 1 divides
     * by zero (producing NaN), and a count below 2 is degenerate for visualization.
     *
     * @param name Parameter name, used in the error message
     * @param value Requested line count
     * @throws std::invalid_argument if value < 2
     */
    static void validateLineCount(const char* name, int value);

    /**
     * @brief Generate a polar grid
     */
    void generatePolarGrid();

    /**
     * @brief Generate a Cartesian grid
     */
    void generateCartesianGrid();

    /**
     * @brief Generate a parametric grid
     */
    void generateParametricGrid();

public:
    /**
     * @brief Construct a new Grid object with default parameters
     *
     * @param gridType Type of grid to create (default: POLAR)
     */
    Grid(GridType gridType = GridType::POLAR);

    /**
     * @brief Construct a new Polar Grid
     *
     * @param numRadial Number of radial lines
     * @param numAngular Number of angular (concentric) lines
     * @param rMinimum Minimum radius (default: 0.0)
     * @param rMaximum Maximum radius (default: 1.0)
     * @param aMinimum Minimum angle in radians (default: 0.0)
     * @param aMaximum Maximum angle in radians (default: 2π)
     * @param includeOriginPoint Whether to include the origin (default: true)
     */
    static Grid createPolarGrid(
        int numRadial,
        int numAngular,
        double rMinimum = 0.0,
        double rMaximum = 1.0,
        double aMinimum = 0.0,
        double aMaximum = 2.0 * M_PI,
        bool includeOriginPoint = true
    );

    /**
     * @brief Construct a new Cartesian Grid
     *
     * @param numHorizontal Number of horizontal lines
     * @param numVertical Number of vertical lines
     * @param xMinimum Minimum x-coordinate (default: -1.0)
     * @param xMaximum Maximum x-coordinate (default: 1.0)
     * @param yMinimum Minimum y-coordinate (default: -1.0)
     * @param yMaximum Maximum y-coordinate (default: 1.0)
     */
    static Grid createCartesianGrid(
        int numHorizontal,
        int numVertical,
        double xMinimum = -1.0,
        double xMaximum = 1.0,
        double yMinimum = -1.0,
        double yMaximum = 1.0
    );

    /**
     * @brief Construct a new Parametric Grid
     *
     * @param paramFunc Function mapping (u,v) to complex points
     * @param numUlines Number of constant-u lines
     * @param numVlines Number of constant-v lines
     * @param uMinimum Minimum u-parameter (default: 0.0)
     * @param uMaximum Maximum u-parameter (default: 1.0)
     * @param vMinimum Minimum v-parameter (default: 0.0)
     * @param vMaximum Maximum v-parameter (default: 1.0)
     */
    static Grid createParametricGrid(
        std::function<Complex(double, double)> paramFunc,
        int numUlines,
        int numVlines,
        double uMinimum = 0.0,
        double uMaximum = 1.0,
        double vMinimum = 0.0,
        double vMaximum = 1.0
    );

    /**
     * @brief Apply a transformation to the grid points
     *
     * @param transform Function mapping complex points to complex points
     * @return Grid New grid with transformed points
     */
    Grid transform(std::function<Complex(const Complex&)> transform) const;

    /**
     * @brief Get all points in the grid
     *
     * @return const std::vector<Complex>& Vector of all grid points
     */
    const std::vector<Complex>& getPoints() const
    {
        return points;
    }

    /**
     * @brief Get line connectivity information
     *
     * Each line is represented as a vector of indices into the points vector
     *
     * @return const std::vector<std::vector<size_t>>& Vector of lines
     */
    const std::vector<std::vector<size_t>>& getLines() const
    {
        return lines;
    }

    /**
     * @brief Get the type of this grid
     *
     * @return GridType Type of the grid
     */
    GridType getType() const
    {
        return type;
    }

    /**
     * @brief Regenerate the grid with current parameters
     */
    void regenerate();

    /**
     * @brief Get the grid parameters as a string map
     *
     * @return std::unordered_map<std::string, double> Map of parameter names to values
     */
    std::unordered_map<std::string, double> getParameters() const;

    /**
     * @brief Set a grid parameter by name
     *
     * @param name Parameter name
     * @param value Parameter value
     * @return bool True if parameter was set, false if name was not recognized
     */
    bool setParameter(const std::string& name, double value);

    /**
     * @brief Add a point to the grid
     *
     * @param point Complex point to add
     * @return size_t Index of the added point
     */
    size_t addPoint(const Complex& point);

    /**
     * @brief Add a line to the grid connecting specified points
     *
     * @param pointIndices Indices of points forming the line
     */
    void addLine(const std::vector<size_t>& pointIndices);

    /**
     * @brief Clear all points and lines from the grid
     */
    void clear();

    /**
     * @brief Get the number of points in the grid
     *
     * @return size_t Number of points
     */
    size_t getPointCount() const
    {
        return points.size();
    }

    /**
     * @brief Get the number of lines in the grid
     *
     * @return size_t Number of lines
     */
    size_t getLineCount() const
    {
        return lines.size();
    }

    /**
     * @brief Create a uniform mesh over the grid
     *
     * Creates a mesh with triangular or quadrilateral elements
     *
     * @param triangulate If true, creates triangular elements, otherwise quadrilateral
     * @return std::vector<std::vector<size_t>> Vector of mesh elements (each element is a vector of point indices)
     */
    std::vector<std::vector<size_t>> createMesh(bool triangulate = true) const;

    /**
     * @brief Create a union of two grids
     *
     * @param other Grid to merge with this one
     * @return Grid Combined grid
     */
    Grid mergeWith(const Grid& other) const;
};

#endif // GRID_HPP
