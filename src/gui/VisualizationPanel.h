#pragma once

#include "../core/Types.h"
#include <memory>
#include <string>
#include <vector>

class ConformalMap;
class Domain;

struct BoundaryCurve
{
    std::string label;
    std::vector<double> x, y;
};

struct GridLine
{
    std::vector<double> x, y;
};

struct BoundingBox
{
    double xMin, xMax, yMin, yMax;
};

/**
 * @brief Panel for visualizing conformal mappings with dual-domain display
 *
 * Provides side-by-side visualization of source and target domains with
 * conformal grid overlays. Uses ImPlot for rendering.
 */
class VisualizationPanel
{
private:
    bool m_showGrid;
    bool m_showSourceDomain;
    bool m_showTargetDomain;
    int m_gridDensity;

    std::vector<BoundaryCurve> m_sourceBoundaries;
    std::vector<BoundaryCurve> m_targetBoundaries;

    std::vector<GridLine> m_sourceGridLines;
    std::vector<GridLine> m_targetGridLines;

    std::shared_ptr<ConformalMap> mp_currentMap;
    bool m_mapComputed;

public:
    VisualizationPanel();
    ~VisualizationPanel();

    bool initialize();
    void render();
    void shutdown();

    /**
     * @brief Update the visualization with a new conformal map
     * @param map Conformal map to visualize
     * @param isComputed Whether the map has been successfully computed. When false,
     *                    target-domain point mapping is skipped (source/target boundary
     *                    geometry is still shown, but no image grid is evaluated) since
     *                    ConformalMap::map() would throw for an uncomputed map.
     */
    void updateMap(std::shared_ptr<ConformalMap> map, bool isComputed = false);

    /**
     * @brief Set grid visibility
     * @param visible Whether to show conformal grid
     */
    void setGridVisible(bool visible) { m_showGrid = visible; }

    /**
     * @brief Set grid density
     * @param density Number of grid lines (must be positive)
     */
    void setGridDensity(int density);

    /**
     * @brief Get current grid visibility
     * @return true if grid is visible
     */
    bool isGridVisible() const { return m_showGrid; }

private:
    void generateBoundariesForDomain(std::shared_ptr<Domain> domain, std::vector<BoundaryCurve>& out);
    void generateSourceGrid();
    void generateTargetGrid();

    void renderSourceDomain();
    void renderTargetDomain();

    BoundingBox computeBounds(const std::vector<BoundaryCurve>& boundaries) const;
    std::string determinePlotTitle(std::shared_ptr<Domain> domain, bool isSource) const;

    void clearGridData();
    void plotGridLines(const std::vector<GridLine>& gridLines);
};
