#include "VisualizationPanel.h"
#include "../core/ConformalMap.h"
#include "../numerics/Grid.h"
#include "imgui.h"
#include "implot.h"
#include <cmath>
#include <spdlog/spdlog.h>

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
    generateSourceBoundary();
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
    
    if (ImPlot::BeginPlot("Source Domain (Unit Circle)", plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");
        ImPlot::SetupAxisLimits(ImAxis_X1, -1.5, 1.5, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -1.5, 1.5, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxesLimits(-1.5, 1.5, -1.5, 1.5, ImGuiCond_FirstUseEver);
        
        // Plot source domain boundary (unit circle)
        if (!m_canonicalBoundaryX.empty() && !m_canonicalBoundaryY.empty())
        {
            ImPlot::PlotLine("Unit Circle", 
                           m_canonicalBoundaryX.data(), 
                           m_canonicalBoundaryY.data(), 
                           static_cast<int>(m_canonicalBoundaryX.size()));
        }
        
        // Plot conformal grid if enabled
        if (m_showGrid && !m_canonicalGridX.empty() && !m_canonicalGridY.empty())
        {
            // Plot grid lines in source domain (unit circle)
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            
            // Plot canonical grid (circles and radials in unit disk)
            ImPlot::PlotLine("Grid", 
                           m_canonicalGridX.data(), 
                           m_canonicalGridY.data(), 
                           static_cast<int>(m_canonicalGridX.size()));
            
            ImPlot::PopStyleColor();
        }
        
        ImPlot::EndPlot();
    }
}

void VisualizationPanel::renderTargetDomain()
{
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = (ImGui::GetContentRegionAvail().x);
    
    if (ImPlot::BeginPlot("Target Domain (Starlike)", plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");
        ImPlot::SetupAxisLimits(ImAxis_X1, -3.0, 3.0, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -2.0, 2.0, ImGuiCond_FirstUseEver);
        ImPlot::SetupAxesLimits(-3.0, 3.0, -2.0, 2.0, ImGuiCond_FirstUseEver);
        
        // Plot target domain boundary (starlike)
        if (!m_targetBoundaryX.empty() && !m_targetBoundaryY.empty())
        {
            ImPlot::PlotLine("Starlike Boundary", 
                           m_targetBoundaryX.data(), 
                           m_targetBoundaryY.data(), 
                           static_cast<int>(m_targetBoundaryX.size()));
        }
        
        // Plot mapped conformal grid if enabled and available
        if (m_showGrid && mp_currentMap && !m_targetGridX.empty() && !m_targetGridY.empty())
        {
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            
            // Plot the mapped grid points as scatter plot
            ImPlot::PlotScatter("Mapped Grid", 
                              m_targetGridX.data(), 
                              m_targetGridY.data(), 
                              static_cast<int>(m_targetGridX.size()));
            
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
    // at grid points in the source domain (ellipse)

    auto sourceDomain = mp_currentMap->getSourceDomain();
    auto starlikeDomain = std::dynamic_pointer_cast<StarlikeDomain>(sourceDomain);
    if (!starlikeDomain)
    {
        return;
    }

    Complex center = starlikeDomain->getCenter();
    int pointsPerLine = 50;

    int totalPoints = 0;
    int failedPoints = 0;
    bool firstFailureLogged = false;

    // Map radial lines from ellipse center to boundary
    for (int i = 0; i < m_gridDensity; ++i)
    {
        double angle = 2.0 * M_PI * i / m_gridDensity;
        double maxRadius = starlikeDomain->getRadius(angle);

        for (int j = 1; j <= pointsPerLine; ++j) // Skip center point
        {
            ++totalPoints;
            double r = static_cast<double>(j) / pointsPerLine * maxRadius;
            Complex z = center + Complex(r * cos(angle), r * sin(angle));

            try
            {
                Complex w = mp_currentMap->map(z);
                m_targetGridX.push_back(w.real());
                m_targetGridY.push_back(w.imag());
            }
            catch (const std::invalid_argument& e)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::warn("Grid radial map evaluation configuration error: {}", e.what());
                    firstFailureLogged = true;
                }
            }
            catch (const std::runtime_error& e)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::debug("Grid radial map evaluation failed at z=({}, {}): {}",
                                  z.real(), z.imag(), e.what());
                    firstFailureLogged = true;
                }
            }
            catch (...)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::warn("Grid radial map evaluation: unknown error at z=({}, {})",
                                 z.real(), z.imag());
                    firstFailureLogged = true;
                }
            }
        }
    }

    // Map elliptical contour lines (scaled ellipses)
    int numContours = m_gridDensity / 2;
    for (int i = 1; i < numContours; ++i)
    {
        double scale = static_cast<double>(i) / numContours;

        for (int j = 0; j <= 100; ++j)
        {
            ++totalPoints;
            double angle = 2.0 * M_PI * j / 100;
            double radius = starlikeDomain->getRadius(angle) * scale;
            Complex z = center + Complex(radius * cos(angle), radius * sin(angle));

            try
            {
                Complex w = mp_currentMap->map(z);
                m_targetGridX.push_back(w.real());
                m_targetGridY.push_back(w.imag());
            }
            catch (const std::invalid_argument& e)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::warn("Grid contour map evaluation configuration error: {}", e.what());
                    firstFailureLogged = true;
                }
            }
            catch (const std::runtime_error& e)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::debug("Grid contour map evaluation failed at z=({}, {}): {}",
                                  z.real(), z.imag(), e.what());
                    firstFailureLogged = true;
                }
            }
            catch (...)
            {
                ++failedPoints;
                if (!firstFailureLogged)
                {
                    spdlog::warn("Grid contour map evaluation: unknown error at z=({}, {})",
                                 z.real(), z.imag());
                    firstFailureLogged = true;
                }
            }
        }
    }

    if (totalPoints > 0 && failedPoints > 0)
    {
        double failureRate = static_cast<double>(failedPoints) / totalPoints;
        if (failureRate > 0.1)
        {
            spdlog::warn("Grid generation: {}/{} points ({:.1f}%) failed map evaluation",
                         failedPoints, totalPoints, failureRate * 100.0);
        }
        else
        {
            spdlog::debug("Grid generation: {}/{} points failed map evaluation",
                          failedPoints, totalPoints);
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

void VisualizationPanel::generateSourceBoundary()
{
    if (!mp_currentMap)
    {
        return;
    }
    
    m_sourceBoundaryX.clear();
    m_sourceBoundaryY.clear();
    
    // Generate source domain boundary from source domain
    auto sourceDomain = mp_currentMap->getSourceDomain();
    if (!sourceDomain)
    {
        return;
    }
    
    // Sample boundary points using domain's getRadius method
    auto starlikeDomain = std::dynamic_pointer_cast<StarlikeDomain>(sourceDomain);
    if (!starlikeDomain)
    {
        return;
    }
    
    Complex center = starlikeDomain->getCenter();
    int numPoints = 200;
    for (int i = 0; i <= numPoints; ++i)
    {
        double angle = 2.0 * M_PI * i / numPoints;
        double radius = starlikeDomain->getRadius(angle);
        
        Complex boundary_point = center + Complex(radius, 0.0) * Complex(cos(angle), sin(angle));
        m_sourceBoundaryX.push_back(boundary_point.real());
        m_sourceBoundaryY.push_back(boundary_point.imag());
    }
}

void VisualizationPanel::generateTargetBoundary()
{
    m_targetBoundaryX.clear();
    m_targetBoundaryY.clear();
    
    // For Theodorsen method, target is always unit circle
    int numPoints = 200;
    for (int i = 0; i <= numPoints; ++i)
    {
        double angle = 2.0 * M_PI * i / numPoints;
        m_targetBoundaryX.push_back(cos(angle));
        m_targetBoundaryY.push_back(sin(angle));
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
