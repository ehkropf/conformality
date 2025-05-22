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

#include <gtest/gtest.h>
#include "../src/Grid.h"
#include <cmath>

// Helper function to compare complex points
bool pointsEqual(const ComplexDouble& a, const ComplexDouble& b, double tolerance = 1e-10)
{
    return std::abs((a - b).getValue()) < tolerance;
}

TEST(GridTest, DefaultConstructionPolar)
{
    Grid grid;

    // Default grid should be polar
    EXPECT_EQ(GridType::POLAR, grid.getType());

    // Check number of lines and points
    EXPECT_TRUE(grid.getPointCount() > 0);
    EXPECT_TRUE(grid.getLineCount() > 0);

    // Default parameters for a polar grid
    auto params = grid.getParameters();
    EXPECT_DOUBLE_EQ(0.0, params["rMin"]);
    EXPECT_DOUBLE_EQ(1.0, params["rMax"]);
    EXPECT_DOUBLE_EQ(0.0, params["angleMin"]);
    EXPECT_DOUBLE_EQ(2.0 * M_PI, params["angleMax"]);
    EXPECT_DOUBLE_EQ(12, params["numRadialLines"]);
    EXPECT_DOUBLE_EQ(5, params["numAngularLines"]);
    EXPECT_DOUBLE_EQ(1.0, params["includeOrigin"]);
}

TEST(GridTest, CreatePolarGrid)
{
    // Create a custom polar grid
    Grid grid = Grid::createPolarGrid(8, 4, 0.5, 2.0, 0.0, M_PI, false);

    EXPECT_EQ(GridType::POLAR, grid.getType());

    // Check parameters
    auto params = grid.getParameters();
    EXPECT_DOUBLE_EQ(0.5, params["rMin"]);
    EXPECT_DOUBLE_EQ(2.0, params["rMax"]);
    EXPECT_DOUBLE_EQ(0.0, params["angleMin"]);
    EXPECT_DOUBLE_EQ(M_PI, params["angleMax"]);
    EXPECT_DOUBLE_EQ(8, params["numRadialLines"]);
    EXPECT_DOUBLE_EQ(4, params["numAngularLines"]);
    EXPECT_DOUBLE_EQ(0.0, params["includeOrigin"]);

    // Check the number of lines
    // Should have 8 radial lines + 4 angular lines = 12 total
    EXPECT_EQ(12, grid.getLineCount());

    // Check if points on the boundary have the correct radius
    const auto& points = grid.getPoints();
    const auto& lines = grid.getLines();

    // Find radial lines (they will have fewer points than angular lines)
    for (const auto& line : lines)
    {
        if (line.size() == 4) // A radial line with 4 points
        {
            // Check the last point (should be at r = 2.0)
            size_t lastIdx = line.back();
            ComplexDouble lastPoint = points[lastIdx];
            EXPECT_NEAR(2.0, lastPoint.abs(), 1e-10);
        }
    }
}

TEST(GridTest, CreateCartesianGrid)
{
    // Create a custom Cartesian grid
    Grid grid = Grid::createCartesianGrid(5, 7, -2.0, 2.0, -1.0, 1.0);

    EXPECT_EQ(GridType::CARTESIAN, grid.getType());

    // Check parameters
    auto params = grid.getParameters();
    EXPECT_DOUBLE_EQ(-2.0, params["xMin"]);
    EXPECT_DOUBLE_EQ(2.0, params["xMax"]);
    EXPECT_DOUBLE_EQ(-1.0, params["yMin"]);
    EXPECT_DOUBLE_EQ(1.0, params["yMax"]);
    EXPECT_DOUBLE_EQ(5, params["numHorizontalLines"]);
    EXPECT_DOUBLE_EQ(7, params["numVerticalLines"]);

    // Check the number of lines
    // Should have 5 horizontal + 7 vertical = 12 total
    EXPECT_EQ(12, grid.getLineCount());

    // Check the number of points
    // Should have 5 * 7 = 35 points
    EXPECT_EQ(35, grid.getPointCount());

    // Check some specific points
    const auto& points = grid.getPoints();

    // Bottom-left corner should be at (-2.0, -1.0)
    EXPECT_NEAR(-2.0, points[0].real(), 1e-10);
    EXPECT_NEAR(-1.0, points[0].imag(), 1e-10);

    // Top-right corner should be at (2.0, 1.0)
    EXPECT_NEAR(2.0, points[34].real(), 1e-10);
    EXPECT_NEAR(1.0, points[34].imag(), 1e-10);
}

TEST(GridTest, CreateParametricGrid)
{
    // Create a parametric grid representing a circle
    auto circleFunc = [](double u, double v) -> ComplexDouble
    {
        // u is angle, v is radius
        return ComplexDouble::fromPolar(v, u);
    };

    Grid grid = Grid::createParametricGrid(circleFunc, 8, 4, 0.0, 2.0 * M_PI, 0.5, 2.0);

    EXPECT_EQ(GridType::PARAMETRIC, grid.getType());

    // Check parameters
    auto params = grid.getParameters();
    EXPECT_DOUBLE_EQ(0.0, params["uMin"]);
    EXPECT_DOUBLE_EQ(2.0 * M_PI, params["uMax"]);
    EXPECT_DOUBLE_EQ(0.5, params["vMin"]);
    EXPECT_DOUBLE_EQ(2.0, params["vMax"]);
    EXPECT_DOUBLE_EQ(8, params["numULines"]);
    EXPECT_DOUBLE_EQ(4, params["numVLines"]);

    // Check the number of lines
    // Should have 8 u-constant lines + 4 v-constant lines = 12 total
    EXPECT_EQ(12, grid.getLineCount());

    // Check the number of points
    // Should have 8 * 4 = 32 points
    EXPECT_EQ(32, grid.getPointCount());

    // Check some specific points
    const auto& points = grid.getPoints();

    // First point should be at angle=0, radius=0.5
    EXPECT_NEAR(0.5, points[0].real(), 1e-10);
    EXPECT_NEAR(0.0, points[0].imag(), 1e-10);

    // Points are stored as [u0,v0], [u0,v1], [u0,v2], [u0,v3], [u1,v0], ...
    // So points with maximum radius (v=2.0) are at indices: 3, 7, 11, 15, 19, 23, 27, 31
    std::vector<int> maxRadiusIndices = {3, 7, 11, 15, 19, 23, 27, 31};
    for (int idx : maxRadiusIndices)
    {
        EXPECT_NEAR(2.0, points[idx].abs(), 1e-10);
    }
}

TEST(GridTest, TransformGrid)
{
    // Create a simple cartesian grid
    Grid grid = Grid::createCartesianGrid(3, 3, -1.0, 1.0, -1.0, 1.0);

    // Define a conformal transformation (z -> z^2)
    auto squareMap = [](const ComplexDouble& z) -> ComplexDouble
    {
        return z * z;
    };

    // Apply the transformation
    Grid transformed = grid.transform(squareMap);

    EXPECT_EQ(grid.getPointCount(), transformed.getPointCount());
    EXPECT_EQ(grid.getLineCount(), transformed.getLineCount());

    // Check that the transformation was applied correctly
    const auto& originalPoints = grid.getPoints();
    const auto& transformedPoints = transformed.getPoints();

    for (size_t i = 0; i < originalPoints.size(); ++i)
    {
        ComplexDouble expected = squareMap(originalPoints[i]);
        EXPECT_TRUE(pointsEqual(expected, transformedPoints[i]));
    }

    // The point at (1,1) should be mapped to (0,2)
    for (size_t i = 0; i < originalPoints.size(); ++i)
    {
        if (std::abs(originalPoints[i].real() - 1.0) < 1e-10 &&
            std::abs(originalPoints[i].imag() - 1.0) < 1e-10)
        {
            EXPECT_NEAR(0.0, transformedPoints[i].real(), 1e-10);
            EXPECT_NEAR(2.0, transformedPoints[i].imag(), 1e-10);
            break;
        }
    }
}

TEST(GridTest, ParameterGetSet)
{
    Grid grid = Grid::createPolarGrid(12, 5);

    // Check initial values
    EXPECT_EQ(12, static_cast<int>(grid.getParameters()["numRadialLines"]));
    EXPECT_EQ(5, static_cast<int>(grid.getParameters()["numAngularLines"]));

    // Set new values
    EXPECT_TRUE(grid.setParameter("numRadialLines", 8.0));
    EXPECT_TRUE(grid.setParameter("numAngularLines", 3.0));

    // Verify the values were changed
    EXPECT_EQ(8, static_cast<int>(grid.getParameters()["numRadialLines"]));
    EXPECT_EQ(3, static_cast<int>(grid.getParameters()["numAngularLines"]));

    // Try setting an invalid parameter
    EXPECT_FALSE(grid.setParameter("invalidParameter", 10.0));

    // Change grid type
    EXPECT_TRUE(grid.setParameter("type", static_cast<double>(GridType::CARTESIAN)));
    EXPECT_EQ(GridType::CARTESIAN, grid.getType());

    // After changing type, should be able to set Cartesian parameters
    EXPECT_TRUE(grid.setParameter("numHorizontalLines", 7.0));
    EXPECT_TRUE(grid.setParameter("numVerticalLines", 9.0));

    // Verify the values were changed
    EXPECT_EQ(7, static_cast<int>(grid.getParameters()["numHorizontalLines"]));
    EXPECT_EQ(9, static_cast<int>(grid.getParameters()["numVerticalLines"]));
}

TEST(GridTest, CreateMesh)
{
    // Create a 3x3 cartesian grid
    Grid grid = Grid::createCartesianGrid(3, 3, 0.0, 2.0, 0.0, 2.0);

    // Create triangular mesh
    auto triangles = grid.createMesh(true);

    // For a 3x3 grid, we should have 2x2 cells, each with 2 triangles
    EXPECT_EQ(8, triangles.size());

    // Each triangle should have exactly 3 vertices
    for (const auto& triangle : triangles)
    {
        EXPECT_EQ(3, triangle.size());
    }

    // Create quadrilateral mesh
    auto quads = grid.createMesh(false);

    // For a 3x3 grid, we should have 2x2 = 4 quads
    EXPECT_EQ(4, quads.size());

    // Each quad should have exactly 4 vertices
    for (const auto& quad : quads)
    {
        EXPECT_EQ(4, quad.size());
    }
}

TEST(GridTest, MergeGrids)
{
    // Create two different grids
    Grid polarGrid = Grid::createPolarGrid(8, 3, 0.0, 1.0);
    Grid cartesianGrid = Grid::createCartesianGrid(3, 3, 1.5, 2.5, 1.5, 2.5);

    size_t polarPoints = polarGrid.getPointCount();
    size_t polarLines = polarGrid.getLineCount();
    size_t cartesianPoints = cartesianGrid.getPointCount();
    size_t cartesianLines = cartesianGrid.getLineCount();

    // Merge the grids
    Grid mergedGrid = polarGrid.mergeWith(cartesianGrid);

    // The merged grid should contain all points and lines from both grids
    EXPECT_EQ(polarPoints + cartesianPoints, mergedGrid.getPointCount());
    EXPECT_EQ(polarLines + cartesianLines, mergedGrid.getLineCount());

    // Check that both original point sets are present in the merged grid
    const auto& polarPoints_orig = polarGrid.getPoints();
    const auto& cartesianPoints_orig = cartesianGrid.getPoints();
    const auto& mergedPoints = mergedGrid.getPoints();

    // First points should match polar grid points
    for (size_t i = 0; i < polarPoints; ++i)
    {
        EXPECT_TRUE(pointsEqual(polarPoints_orig[i], mergedPoints[i]));
    }

    // Later points should match cartesian grid points
    for (size_t i = 0; i < cartesianPoints; ++i)
    {
        EXPECT_TRUE(pointsEqual(cartesianPoints_orig[i], mergedPoints[i + polarPoints]));
    }
}

TEST(GridTest, AddingPointsAndLines)
{
    Grid grid;
    grid.clear();

    // Add individual points
    size_t idx1 = grid.addPoint(ComplexDouble(0.0, 0.0));
    size_t idx2 = grid.addPoint(ComplexDouble(1.0, 0.0));
    size_t idx3 = grid.addPoint(ComplexDouble(0.0, 1.0));

    EXPECT_EQ(0, idx1);
    EXPECT_EQ(1, idx2);
    EXPECT_EQ(2, idx3);
    EXPECT_EQ(3, grid.getPointCount());

    // Add a line connecting the points
    grid.addLine({idx1, idx2, idx3, idx1});

    EXPECT_EQ(1, grid.getLineCount());

    // Check the line's point indices
    const auto& lines = grid.getLines();
    EXPECT_EQ(4, lines[0].size());
    EXPECT_EQ(idx1, lines[0][0]);
    EXPECT_EQ(idx2, lines[0][1]);
    EXPECT_EQ(idx3, lines[0][2]);
    EXPECT_EQ(idx1, lines[0][3]);

    // Adding a line with an invalid index should throw
    EXPECT_THROW(grid.addLine({0, 1, 99}), std::out_of_range);
}

TEST(GridTest, Regeneration)
{
    Grid grid = Grid::createPolarGrid(8, 4);

    size_t initialPoints = grid.getPointCount();
    size_t initialLines = grid.getLineCount();

    // Modify parameters and regenerate
    grid.setParameter("numRadialLines", 12.0);
    grid.setParameter("numAngularLines", 6.0);
    grid.regenerate();

    // Point and line counts should change
    EXPECT_NE(initialPoints, grid.getPointCount());
    EXPECT_NE(initialLines, grid.getLineCount());

    // Parameters should be updated
    EXPECT_EQ(12, static_cast<int>(grid.getParameters()["numRadialLines"]));
    EXPECT_EQ(6, static_cast<int>(grid.getParameters()["numAngularLines"]));

    // Number of lines should be 12 radial + (6-1) angular = 17 total
    // (One angular line is skipped because it would be degenerate at r=0)
    EXPECT_EQ(17, grid.getLineCount());
}
