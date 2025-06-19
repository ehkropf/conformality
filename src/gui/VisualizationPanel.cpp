#include "VisualizationPanel.h"
#include "../ConformalMap.h"
#include "../Grid.h"
#include "imgui.h"
#include "implot.h"
#include <cmath>

VisualizationPanel::VisualizationPanel()
    : m_showGrid{true}
    , m_showCanonicalDomain{true}
    , m_showTargetDomain{true}
    , m_gridDensity{8}
    , mp_currentMap{nullptr}
{
}

VisualizationPanel::~VisualizationPanel()
{
    shutdown();
}

bool VisualizationPanel::initialize()
{
    // Generate default unit circle boundary
    generateCanonicalBoundary();
    generateCanonicalGrid();
    
    return true;
}

void VisualizationPanel::render()
{
    // Split the visualization area into two plots side by side
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
    
    // Canonical domain plot (left side)
    if (m_showCanonicalDomain)
    {
        renderCanonicalDomain();
    }
    
    ImGui::SameLine();
    
    // Target domain plot (right side)
    if (m_showTargetDomain)
    {
        renderTargetDomain();
    }
}

void VisualizationPanel::shutdown()
{
    clearGridData();
    mp_currentMap.reset();
}

void VisualizationPanel::updateMap(std::shared_ptr<ConformalMap> map)
{
    mp_currentMap = map;
    
    // Regenerate visualization data
    generateTargetBoundary();
    if (m_showGrid)
    {
        generateTargetGrid();
    }
}

void VisualizationPanel::setGridDensity(int density)
{
    if (density > 0 && density != m_gridDensity)
    {
        m_gridDensity = density;
        generateCanonicalGrid();
        if (mp_currentMap)
        {
            generateTargetGrid();
        }
    }
}

void VisualizationPanel::renderCanonicalDomain()
{
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
    
    if (ImPlot::BeginPlot("Canonical Domain (Unit Circle)", plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.5, 1.5, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.5, 1.5, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxesLimits(-1.5, 1.5, -1.5, 1.5, ImGuiCond_FirstUseEver);
        
        // Plot unit circle boundary
        if (!m_canonicalBoundaryX.empty() && !m_canonicalBoundaryY.empty())
        {
            ImPlot::PlotLine("Unit Circle", 
                           m_canonicalBoundaryX.data(), 
                           m_canonicalBoundaryY.data(), 
                           static_cast<int>(m_canonicalBoundaryX.size()));
        }
        
        // Plot conformal grid if enabled
        if (m_showGrid)
        {
            // Plot grid lines
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            
            // Radial grid lines (from center to boundary)
            int pointsPerLine = 50;
            for (int i = 0; i < m_gridDensity; ++i)
            {
                double angle = 2.0 * M_PI * i / m_gridDensity;
                std::vector<double> radialX, radialY;
                
                for (int j = 0; j <= pointsPerLine; ++j)
                {
                    double r = static_cast<double>(j) / pointsPerLine;
                    radialX.push_back(r * cos(angle));
                    radialY.push_back(r * sin(angle));
                }
                
                ImPlot::PlotLine(("Radial" + std::to_string(i)).c_str(), 
                               radialX.data(), radialY.data(), radialX.size());
            }
            
            // Circular grid lines (concentric circles)
            int numCircles = m_gridDensity / 2;
            for (int i = 1; i < numCircles; ++i)
            {
                double radius = static_cast<double>(i) / numCircles;
                std::vector<double> circleX, circleY;
                
                for (int j = 0; j <= 100; ++j)
                {
                    double angle = 2.0 * M_PI * j / 100;
                    circleX.push_back(radius * cos(angle));
                    circleY.push_back(radius * sin(angle));
                }
                
                ImPlot::PlotLine(("Circle" + std::to_string(i)).c_str(), 
                               circleX.data(), circleY.data(), circleX.size());
            }
            
            ImPlot::PopStyleColor();
        }
        
        ImPlot::EndPlot();
    }
}

void VisualizationPanel::renderTargetDomain()
{
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = (ImGui::GetContentRegionAvail().x);
    
    if (ImPlot::BeginPlot("Target Domain", plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");
        ImPlot::SetupAxisLimits(ImAxis_X1, -3.0, 3.0, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -2.0, 2.0, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxesLimits(-3.0, 3.0, -2.0, 2.0, ImGuiCond_FirstUseEver);
        
        // Plot target domain boundary
        if (!m_targetBoundaryX.empty() && !m_targetBoundaryY.empty())
        {
            ImPlot::PlotLine("Target Boundary", 
                           m_targetBoundaryX.data(), 
                           m_targetBoundaryY.data(), 
                           static_cast<int>(m_targetBoundaryX.size()));
        }
        
        // Plot mapped conformal grid if enabled and available
        if (m_showGrid && mp_currentMap && !m_targetGridX.empty() && !m_targetGridY.empty())
        {
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            
            // This would plot the mapped grid lines
            // Implementation depends on how the grid data is structured
            // For now, we'll show placeholder text
            
            ImPlot::PopStyleColor();
        }
        
        // Show message if no map is loaded
        if (!mp_currentMap)
        {
            ImPlot::PlotText("No mapping computed", 0.0, 0.0);
        }
        
        ImPlot::EndPlot();
    }
}

void VisualizationPanel::generateCanonicalGrid()
{
    clearGridData();
    
    // Generate grid data for the canonical domain (unit circle)
    // Pre-generate data to avoid doing it every frame in render
    
    // We don't actually need to store the grid data since it's generated
    // algorithmically in renderCanonicalDomain(). The grid visibility
    // is controlled by m_showGrid flag during rendering.
}

void VisualizationPanel::generateTargetGrid()
{
    if (!mp_currentMap)
    {
        return;
    }
    
    m_targetGridX.clear();
    m_targetGridY.clear();
    
    // Generate mapped grid points by evaluating the conformal map
    // at grid points in the canonical domain
    
    int pointsPerLine = 50;
    
    // Map radial lines
    for (int i = 0; i < m_gridDensity; ++i)
    {
        double angle = 2.0 * M_PI * i / m_gridDensity;
        
        for (int j = 1; j <= pointsPerLine; ++j) // Skip center point
        {
            double r = static_cast<double>(j) / pointsPerLine;
            Complex z(r * cos(angle), r * sin(angle));
            
            try
            {
                Complex w = mp_currentMap->map(z);
                m_targetGridX.push_back(w.real());
                m_targetGridY.push_back(w.imag());
            }
            catch (...)
            {
                // Skip points where evaluation fails
            }
        }
    }
    
    // Map circular lines
    int numCircles = m_gridDensity / 2;
    for (int i = 1; i < numCircles; ++i)
    {
        double radius = static_cast<double>(i) / numCircles;
        
        for (int j = 0; j <= 100; ++j)
        {
            double angle = 2.0 * M_PI * j / 100;
            Complex z(radius * cos(angle), radius * sin(angle));
            
            try
            {
                Complex w = mp_currentMap->map(z);
                m_targetGridX.push_back(w.real());
                m_targetGridY.push_back(w.imag());
            }
            catch (...)
            {
                // Skip points where evaluation fails
            }
        }
    }
}

void VisualizationPanel::generateCanonicalBoundary()
{
    m_canonicalBoundaryX.clear();
    m_canonicalBoundaryY.clear();
    
    // Generate unit circle boundary
    int numPoints = 200;
    for (int i = 0; i <= numPoints; ++i)
    {
        double angle = 2.0 * M_PI * i / numPoints;
        m_canonicalBoundaryX.push_back(cos(angle));
        m_canonicalBoundaryY.push_back(sin(angle));
    }
}

void VisualizationPanel::generateTargetBoundary()
{
    if (!mp_currentMap)
    {
        return;
    }
    
    m_targetBoundaryX.clear();
    m_targetBoundaryY.clear();
    
    // Generate target boundary by mapping unit circle points
    int numPoints = 200;
    for (int i = 0; i <= numPoints; ++i)
    {
        double angle = 2.0 * M_PI * i / numPoints;
        Complex z(cos(angle), sin(angle));
        
        try
        {
            Complex w = mp_currentMap->map(z);
            m_targetBoundaryX.push_back(w.real());
            m_targetBoundaryY.push_back(w.imag());
        }
        catch (...)
        {
            // Use unit circle point if evaluation fails
            m_targetBoundaryX.push_back(cos(angle));
            m_targetBoundaryY.push_back(sin(angle));
        }
    }
}

void VisualizationPanel::clearGridData()
{
    m_canonicalGridX.clear();
    m_canonicalGridY.clear();
    m_targetGridX.clear();
    m_targetGridY.clear();
    
    // Note: Don't clear boundary data here - it should persist
    // Boundary data is managed separately in generate*Boundary() methods
}