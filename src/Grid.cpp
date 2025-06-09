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

#include "Grid.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

Grid::Grid(GridType gridType)
    : type(gridType)
    , xMin(-1.0)
    , xMax(1.0)
    , yMin(-1.0)
    , yMax(1.0)
    , rMin(0.0)
    , rMax(1.0)
    , angleMin(0.0)
    , angleMax(2.0 * M_PI)
    , numRadialLines(12)
    , numAngularLines(5)
    , numHorizontalLines(11)
    , numVerticalLines(11)
    , includeOrigin(true)
    , uMin(0.0)
    , uMax(1.0)
    , vMin(0.0)
    , vMax(1.0)
    , numULines(10)
    , numVLines(10)
{
    // Default parameterization for parametric grids (identity function)
    parameterization = [](double u, double v) -> Complex
    {
        return Complex(u, v);
    };

    regenerate();
}

Grid Grid::createPolarGrid(
    int numRadial,
    int numAngular,
    double rMinimum,
    double rMaximum,
    double aMinimum,
    double aMaximum,
    bool includeOriginPoint)
{
    Grid grid(GridType::POLAR);

    grid.numRadialLines = numRadial;
    grid.numAngularLines = numAngular;
    grid.rMin = rMinimum;
    grid.rMax = rMaximum;
    grid.angleMin = aMinimum;
    grid.angleMax = aMaximum;
    grid.includeOrigin = includeOriginPoint;

    grid.regenerate();
    return grid;
}

Grid Grid::createCartesianGrid(
    int numHorizontal,
    int numVertical,
    double xMinimum,
    double xMaximum,
    double yMinimum,
    double yMaximum)
{
    Grid grid(GridType::CARTESIAN);

    grid.numHorizontalLines = numHorizontal;
    grid.numVerticalLines = numVertical;
    grid.xMin = xMinimum;
    grid.xMax = xMaximum;
    grid.yMin = yMinimum;
    grid.yMax = yMaximum;

    grid.regenerate();
    return grid;
}

Grid Grid::createParametricGrid(
    std::function<Complex(double, double)> paramFunc,
    int numUlines,
    int numVlines,
    double uMinimum,
    double uMaximum,
    double vMinimum,
    double vMaximum)
{
    Grid grid(GridType::PARAMETRIC);

    grid.parameterization = paramFunc;
    grid.numULines = numUlines;
    grid.numVLines = numVlines;
    grid.uMin = uMinimum;
    grid.uMax = uMaximum;
    grid.vMin = vMinimum;
    grid.vMax = vMaximum;

    grid.regenerate();
    return grid;
}

void Grid::generatePolarGrid()
{
    clear();

    // Add origin if requested
    size_t originIndex = 0;
    if (includeOrigin && rMin <= 0.0)
    {
        originIndex = addPoint(Complex(0.0, 0.0));
    }

    // Create radial lines
    for (int i = 0; i < numRadialLines; ++i)
    {
        double angle = angleMin + i * (angleMax - angleMin) / numRadialLines;

        std::vector<size_t> line;
        if (includeOrigin && rMin <= 0.0)
        {
            line.push_back(originIndex);
        }

        // Add points along this radial line
        for (int j = 0; j < numAngularLines; ++j)
        {
            double r = rMin + j * (rMax - rMin) / (numAngularLines - 1);
            if (r <= 0.0 && j == 0)
            {
                // Skip origin if already added
                if (includeOrigin)
                {
                    continue;
                }
                r = std::min(r, 1e-10); // Avoid exact zero to prevent degenerate point
            }

            Complex point = std::polar(r, angle);
            line.push_back(addPoint(point));
        }

        addLine(line);
    }

    // Create angular (concentric) lines
    for (int j = 0; j < numAngularLines; ++j)
    {
        double r = rMin + j * (rMax - rMin) / (numAngularLines - 1);
        if (r <= 0.0 && j == 0)
        {
            // Skip the innermost angular line if it would be degenerate
            continue;
        }

        std::vector<size_t> line;
        std::vector<Complex> circlePoints;

        // Create points along this circle
        for (int i = 0; i <= numRadialLines; ++i)
        {
            double angle = angleMin + i * (angleMax - angleMin) / numRadialLines;
            if (i == numRadialLines && std::abs(angleMax - angleMin - 2.0 * M_PI) < 1e-10)
            {
                // Skip the last point if it overlaps with the first for a full circle
                break;
            }

            Complex point = std::polar(r, angle);
            circlePoints.push_back(point);
        }

        // Find or add these points and create the line
        for (const auto& point : circlePoints)
        {
            // Check if this point already exists in the grid
            auto it = std::find_if(points.begin(), points.end(), [&point](const Complex& p)
            {
                return std::abs(point - p) < 1e-10;
            });

            if (it != points.end())
            {
                line.push_back(std::distance(points.begin(), it));
            }
            else
            {
                line.push_back(addPoint(point));
            }
        }

        addLine(line);
    }
}

void Grid::generateCartesianGrid()
{
    clear();

    // Create a 2D array to store point indices for efficient lookup
    std::vector<std::vector<size_t>> pointIndices(numHorizontalLines,
                                                 std::vector<size_t>(numVerticalLines));

    // Create all grid points and store their indices
    for (int j = 0; j < numHorizontalLines; ++j)
    {
        double y = yMin + j * (yMax - yMin) / (numHorizontalLines - 1);

        for (int i = 0; i < numVerticalLines; ++i)
        {
            double x = xMin + i * (xMax - xMin) / (numVerticalLines - 1);
            pointIndices[j][i] = addPoint(Complex(x, y));
        }
    }

    // Create horizontal lines
    for (int j = 0; j < numHorizontalLines; ++j)
    {
        std::vector<size_t> line;

        for (int i = 0; i < numVerticalLines; ++i)
        {
            line.push_back(pointIndices[j][i]);
        }

        addLine(line);
    }

    // Create vertical lines
    for (int i = 0; i < numVerticalLines; ++i)
    {
        std::vector<size_t> line;

        for (int j = 0; j < numHorizontalLines; ++j)
        {
            line.push_back(pointIndices[j][i]);
        }

        addLine(line);
    }
}

void Grid::generateParametricGrid()
{
    clear();

    // Create matrix of points for efficient indexing
    std::vector<std::vector<size_t>> pointIndices(numULines, std::vector<size_t>(numVLines));

    // Create all grid points
    for (int i = 0; i < numULines; ++i)
    {
        double u = uMin + i * (uMax - uMin) / (numULines - 1);

        for (int j = 0; j < numVLines; ++j)
        {
            double v = vMin + j * (vMax - vMin) / (numVLines - 1);
            Complex point = parameterization(u, v);
            pointIndices[i][j] = addPoint(point);
        }
    }

    // Create constant-u lines
    for (int i = 0; i < numULines; ++i)
    {
        std::vector<size_t> line;

        for (int j = 0; j < numVLines; ++j)
        {
            line.push_back(pointIndices[i][j]);
        }

        addLine(line);
    }

    // Create constant-v lines
    for (int j = 0; j < numVLines; ++j)
    {
        std::vector<size_t> line;

        for (int i = 0; i < numULines; ++i)
        {
            line.push_back(pointIndices[i][j]);
        }

        addLine(line);
    }
}

void Grid::regenerate()
{
    switch (type)
    {
        case GridType::POLAR:
            generatePolarGrid();
            break;

        case GridType::CARTESIAN:
            generateCartesianGrid();
            break;

        case GridType::PARAMETRIC:
            generateParametricGrid();
            break;
    }
}

Grid Grid::transform(std::function<Complex(const Complex&)> transform) const
{
    Grid transformedGrid = *this;
    transformedGrid.clear();

    // Transform all points
    std::vector<size_t> newIndices(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    {
        Complex transformed = transform(points[i]);
        newIndices[i] = transformedGrid.addPoint(transformed);
    }

    // Recreate the lines with the new point indices
    for (const auto& line : lines)
    {
        std::vector<size_t> newLine;
        for (size_t pointIndex : line)
        {
            newLine.push_back(newIndices[pointIndex]);
        }
        transformedGrid.addLine(newLine);
    }

    return transformedGrid;
}

size_t Grid::addPoint(const Complex& point)
{
    points.push_back(point);
    return points.size() - 1;
}

void Grid::addLine(const std::vector<size_t>& pointIndices)
{
    if (pointIndices.empty())
    {
        return;
    }

    // Validate that all indices are within bounds
    for (size_t index : pointIndices)
    {
        if (index >= points.size())
        {
            throw std::out_of_range("Point index out of range in addLine");
        }
    }

    lines.push_back(pointIndices);
}

void Grid::clear()
{
    points.clear();
    lines.clear();
}

std::unordered_map<std::string, double> Grid::getParameters() const
{
    std::unordered_map<std::string, double> params;

    // Common parameters
    params["type"] = static_cast<double>(type);

    // Type-specific parameters
    switch (type)
    {
        case GridType::POLAR:
            params["rMin"] = rMin;
            params["rMax"] = rMax;
            params["angleMin"] = angleMin;
            params["angleMax"] = angleMax;
            params["numRadialLines"] = numRadialLines;
            params["numAngularLines"] = numAngularLines;
            params["includeOrigin"] = includeOrigin ? 1.0 : 0.0;
            break;

        case GridType::CARTESIAN:
            params["xMin"] = xMin;
            params["xMax"] = xMax;
            params["yMin"] = yMin;
            params["yMax"] = yMax;
            params["numHorizontalLines"] = numHorizontalLines;
            params["numVerticalLines"] = numVerticalLines;
            break;

        case GridType::PARAMETRIC:
            params["uMin"] = uMin;
            params["uMax"] = uMax;
            params["vMin"] = vMin;
            params["vMax"] = vMax;
            params["numULines"] = numULines;
            params["numVLines"] = numVLines;
            break;
    }

    return params;
}

bool Grid::setParameter(const std::string& name, double value)
{
    // Common parameters
    if (name == "type")
    {
        int typeValue = static_cast<int>(value);
        if (typeValue >= 0 && typeValue <= 2)
        {
            type = static_cast<GridType>(typeValue);
            return true;
        }
        return false;
    }

    // Type-specific parameters
    switch (type)
    {
        case GridType::POLAR:
            if (name == "rMin") { rMin = value; return true; }
            else if (name == "rMax") { rMax = value; return true; }
            else if (name == "angleMin") { angleMin = value; return true; }
            else if (name == "angleMax") { angleMax = value; return true; }
            else if (name == "numRadialLines") { numRadialLines = static_cast<int>(value); return true; }
            else if (name == "numAngularLines") { numAngularLines = static_cast<int>(value); return true; }
            else if (name == "includeOrigin") { includeOrigin = value != 0.0; return true; }
            break;

        case GridType::CARTESIAN:
            if (name == "xMin") { xMin = value; return true; }
            else if (name == "xMax") { xMax = value; return true; }
            else if (name == "yMin") { yMin = value; return true; }
            else if (name == "yMax") { yMax = value; return true; }
            else if (name == "numHorizontalLines") { numHorizontalLines = static_cast<int>(value); return true; }
            else if (name == "numVerticalLines") { numVerticalLines = static_cast<int>(value); return true; }
            break;

        case GridType::PARAMETRIC:
            if (name == "uMin") { uMin = value; return true; }
            else if (name == "uMax") { uMax = value; return true; }
            else if (name == "vMin") { vMin = value; return true; }
            else if (name == "vMax") { vMax = value; return true; }
            else if (name == "numULines") { numULines = static_cast<int>(value); return true; }
            else if (name == "numVLines") { numVLines = static_cast<int>(value); return true; }
            break;
    }

    return false;
}

std::vector<std::vector<size_t>> Grid::createMesh(bool triangulate) const
{
    std::vector<std::vector<size_t>> elements;

    // Only supports Cartesian and parametric grids for now
    if (type != GridType::CARTESIAN && type != GridType::PARAMETRIC)
    {
        return elements;
    }

    int numRows, numCols;
    if (type == GridType::CARTESIAN)
    {
        numRows = numHorizontalLines;
        numCols = numVerticalLines;
    }
    else // PARAMETRIC
    {
        numRows = numULines;
        numCols = numVLines;
    }

    // Create elements (either triangles or quads)
    for (int i = 0; i < numRows - 1; ++i)
    {
        for (int j = 0; j < numCols - 1; ++j)
        {
            // Compute indices of the four corners of this cell
            size_t bottomLeft = i * numCols + j;
            size_t bottomRight = bottomLeft + 1;
            size_t topLeft = (i + 1) * numCols + j;
            size_t topRight = topLeft + 1;

            if (triangulate)
            {
                // Create two triangles for this cell
                elements.push_back({bottomLeft, bottomRight, topLeft});
                elements.push_back({bottomRight, topRight, topLeft});
            }
            else
            {
                // Create one quad for this cell
                elements.push_back({bottomLeft, bottomRight, topRight, topLeft});
            }
        }
    }

    return elements;
}

Grid Grid::mergeWith(const Grid& other) const
{
    Grid mergedGrid = *this;

    // Get the current number of points in this grid (offset for the other grid)
    size_t offset = points.size();

    // Add all points from the other grid
    for (const auto& point : other.points)
    {
        mergedGrid.addPoint(point);
    }

    // Add all lines from the other grid (adjusting indices)
    for (const auto& line : other.lines)
    {
        std::vector<size_t> newLine;
        for (size_t index : line)
        {
            newLine.push_back(index + offset);
        }
        mergedGrid.addLine(newLine);
    }

    return mergedGrid;
}
