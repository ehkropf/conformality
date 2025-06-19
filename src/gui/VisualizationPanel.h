#pragma once

#include "../core/Types.h"
#include <vector>
#include <memory>

class ConformalMap;
class Grid;

/**
 * @brief Panel for visualizing conformal mappings with dual-domain display
 * 
 * Provides side-by-side visualization of the canonical domain (unit circle)
 * and target domain with conformal grid overlays. Uses ImPlot for rendering.
 */
class VisualizationPanel
{
private:
    bool m_showGrid;
    bool m_showCanonicalDomain;
    bool m_showTargetDomain;
    int m_gridDensity;
    
    // Grid data for plotting
    std::vector<double> m_canonicalGridX;
    std::vector<double> m_canonicalGridY;
    std::vector<double> m_targetGridX;
    std::vector<double> m_targetGridY;
    
    // Domain boundary data
    std::vector<double> m_canonicalBoundaryX;
    std::vector<double> m_canonicalBoundaryY;
    std::vector<double> m_sourceBoundaryX;
    std::vector<double> m_sourceBoundaryY;
    std::vector<double> m_targetBoundaryX;
    std::vector<double> m_targetBoundaryY;
    
    std::shared_ptr<ConformalMap> mp_currentMap;
    
public:
    VisualizationPanel();
    ~VisualizationPanel();
    
    bool initialize();
    void render();
    void shutdown();
    
    /**
     * @brief Update the visualization with a new conformal map
     * @param map Conformal map to visualize
     */
    void updateMap(std::shared_ptr<ConformalMap> map);
    
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
    void renderCanonicalDomain();
    void renderTargetDomain();
    void generateCanonicalGrid();
    void generateTargetGrid();
    void generateCanonicalBoundary();
    void generateSourceBoundary();
    void generateTargetBoundary();
    void clearGridData();
};