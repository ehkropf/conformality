#include "VisualizationPanel.h"
#include "../core/ConformalMap.h"
#include "../domains/Domain.h"
#include "../domains/FornbergCanonicalDomain.h"
#include "imgui.h"
#include "implot.h"
#include <cmath>
#include <limits>
#include <spdlog/spdlog.h>

VisualizationPanel::VisualizationPanel()
    : m_showGrid{true}
    , m_showSourceDomain{true}
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
    // Generate default unit circle as placeholder source boundary
    BoundaryCurve unitCircle;
    unitCircle.label = "Unit Circle";
    int numPoints = 200;
    for (int i = 0; i <= numPoints; ++i)
    {
        double angle = 2.0 * M_PI * i / numPoints;
        unitCircle.x.push_back(cos(angle));
        unitCircle.y.push_back(sin(angle));
    }
    m_sourceBoundaries.clear();
    m_sourceBoundaries.push_back(std::move(unitCircle));

    return true;
}

void VisualizationPanel::render()
{
    if (m_showSourceDomain)
    {
        renderSourceDomain();
    }

    ImGui::SameLine();

    if (m_showTargetDomain)
    {
        renderTargetDomain();
    }
}

void VisualizationPanel::shutdown()
{
    clearGridData();
    m_sourceBoundaries.clear();
    m_targetBoundaries.clear();
    mp_currentMap.reset();
}

void VisualizationPanel::updateMap(std::shared_ptr<ConformalMap> map)
{
    mp_currentMap = map;

    if (!mp_currentMap)
    {
        return;
    }

    // Generate boundaries from source and target domains
    auto sourceDomain = mp_currentMap->getSourceDomain();
    auto targetDomain = mp_currentMap->getTargetDomain();

    if (sourceDomain)
    {
        m_sourceBoundaries.clear();
        generateBoundariesForDomain(sourceDomain, m_sourceBoundaries);
    }

    if (targetDomain)
    {
        m_targetBoundaries.clear();
        generateBoundariesForDomain(targetDomain, m_targetBoundaries);
    }

    // Generate grids
    if (sourceDomain)
    {
        generateSourceGrid();
    }
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
        if (mp_currentMap)
        {
            generateSourceGrid();
            generateTargetGrid();
        }
    }
}

// --- Boundary generation ---

void VisualizationPanel::generateBoundariesForDomain(std::shared_ptr<Domain> domain, std::vector<BoundaryCurve>& out)
{
    constexpr int numPoints = 200;

    // Try multiply-connected first (FornbergCanonicalDomain inherits from MC, not SC)
    auto mcDomain = std::dynamic_pointer_cast<MultiplyConnectedDomain>(domain);
    if (mcDomain)
    {
        const auto& boundaries = mcDomain->getBoundaries();
        for (size_t b = 0; b < boundaries.size(); ++b)
        {
            std::vector<std::vector<Complex>> samples;
            try
            {
                samples = boundaries[b]->sample(numPoints);
            }
            catch (const std::exception& e)
            {
                spdlog::warn("Failed to sample boundary {}: {}", b, e.what());
                continue;
            }

            for (size_t comp = 0; comp < samples.size(); ++comp)
            {
                BoundaryCurve curve;
                if (b == 0)
                {
                    curve.label = "Outer Boundary";
                }
                else
                {
                    curve.label = "Hole " + std::to_string(b);
                }

                for (const auto& z : samples[comp])
                {
                    curve.x.push_back(z.real());
                    curve.y.push_back(z.imag());
                }
                // Close the curve
                if (!samples[comp].empty())
                {
                    curve.x.push_back(samples[comp].front().real());
                    curve.y.push_back(samples[comp].front().imag());
                }
                out.push_back(std::move(curve));
            }
        }
        return;
    }

    // Try simply-connected
    auto scDomain = std::dynamic_pointer_cast<SimplyConnectedDomain>(domain);
    if (scDomain)
    {
        std::vector<std::vector<Complex>> samples;
        try
        {
            samples = scDomain->getBoundary().sample(numPoints);
        }
        catch (const std::exception& e)
        {
            spdlog::warn("Failed to sample boundary: {}", e.what());
            return;
        }

        for (size_t comp = 0; comp < samples.size(); ++comp)
        {
            BoundaryCurve curve;
            curve.label = "Boundary";

            for (const auto& z : samples[comp])
            {
                curve.x.push_back(z.real());
                curve.y.push_back(z.imag());
            }
            // Close the curve
            if (!samples[comp].empty())
            {
                curve.x.push_back(samples[comp].front().real());
                curve.y.push_back(samples[comp].front().imag());
            }
            out.push_back(std::move(curve));
        }
        return;
    }

    spdlog::warn("generateBoundariesForDomain: unrecognized domain type; no boundaries generated");
}

// --- Bounds and titles ---

BoundingBox VisualizationPanel::computeBounds(const std::vector<BoundaryCurve>& boundaries) const
{
    BoundingBox box{
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::lowest()
    };

    bool hasData = false;
    for (const auto& curve : boundaries)
    {
        for (size_t i = 0; i < curve.x.size(); ++i)
        {
            double px = curve.x[i];
            double py = curve.y[i];
            if (std::isnan(px) || std::isnan(py))
            {
                continue;
            }
            if (px < box.xMin) box.xMin = px;
            if (px > box.xMax) box.xMax = px;
            if (py < box.yMin) box.yMin = py;
            if (py > box.yMax) box.yMax = py;
            hasData = true;
        }
    }

    if (!hasData)
    {
        box = {-1.5, 1.5, -1.5, 1.5};
    }

    return box;
}

std::string VisualizationPanel::determinePlotTitle(std::shared_ptr<Domain> domain, bool isSource) const
{
    if (!domain)
    {
        return isSource ? "Source Domain" : "Target Domain";
    }

    // Check most specific types first
    if (std::dynamic_pointer_cast<FornbergCanonicalDomain>(domain))
    {
        return "Canonical Domain";
    }
    if (std::dynamic_pointer_cast<CircularDomain>(domain))
    {
        return isSource ? "Source Domain (Circle)" : "Target Domain (Circle)";
    }
    if (std::dynamic_pointer_cast<StarlikeDomain>(domain))
    {
        return isSource ? "Source Domain (Starlike)" : "Target Domain (Starlike)";
    }
    if (std::dynamic_pointer_cast<MultiplyConnectedDomain>(domain))
    {
        return isSource ? "Source Domain (MC)" : "Target Domain (MC)";
    }

    return isSource ? "Source Domain" : "Target Domain";
}

// --- Grid generation ---

void VisualizationPanel::generateSourceGrid()
{
    m_sourceGridLines.clear();

    if (!mp_currentMap)
    {
        return;
    }

    auto sourceDomain = mp_currentMap->getSourceDomain();
    if (!sourceDomain)
    {
        return;
    }

    // For StarlikeDomain, generate polar grid (radials + contours) to reflect star-shaped geometry
    auto starlikeDomain = std::dynamic_pointer_cast<StarlikeDomain>(sourceDomain);
    if (starlikeDomain)
    {
        Complex center = starlikeDomain->getCenter();
        int pointsPerLine = 100;

        // Radial lines
        for (int i = 0; i < m_gridDensity; ++i)
        {
            double angle = 2.0 * M_PI * i / m_gridDensity;
            double maxRadius = starlikeDomain->getRadius(angle);

            GridLine line;
            for (int j = 0; j <= pointsPerLine; ++j)
            {
                double r = static_cast<double>(j) / pointsPerLine * maxRadius;
                Complex z = center + Complex(r * cos(angle), r * sin(angle));
                line.x.push_back(z.real());
                line.y.push_back(z.imag());
            }
            m_sourceGridLines.push_back(std::move(line));
        }

        // Contour lines (scaled boundary curves)
        int numContours = m_gridDensity / 2;
        for (int i = 1; i < numContours; ++i)
        {
            double scale = static_cast<double>(i) / numContours;
            GridLine line;
            for (int j = 0; j <= pointsPerLine; ++j)
            {
                double angle = 2.0 * M_PI * j / pointsPerLine;
                double radius = starlikeDomain->getRadius(angle) * scale;
                Complex z = center + Complex(radius * cos(angle), radius * sin(angle));
                line.x.push_back(z.real());
                line.y.push_back(z.imag());
            }
            m_sourceGridLines.push_back(std::move(line));
        }
        return;
    }

    // For FornbergCanonicalDomain (unit disk with circular holes), use a polar grid
    // with the same geometry as the MATLAB reference (bdd_plot2.m): radial lines +
    // concentric circles. Hole clipping (NaN insertion) is added here since C++ generates
    // the grid directly rather than inverting from the slit domain.
    auto canonicalDomain = std::dynamic_pointer_cast<FornbergCanonicalDomain>(sourceDomain);
    if (canonicalDomain)
    {
        int pointsPerLine = 200;
        int numRadials = 2 * m_gridDensity;
        int numCircles = m_gridDensity;
        const auto& holeCenters = canonicalDomain->getHoleCenters();
        const auto& holeRadii = canonicalDomain->getHoleRadii();

        // Append point to grid line, replacing with NaN if inside any hole
        auto appendPoint = [&holeCenters, &holeRadii](GridLine& line, double x, double y)
        {
            for (size_t h = 0; h < holeCenters.size(); ++h)
            {
                double dx = x - holeCenters[h].real();
                double dy = y - holeCenters[h].imag();
                if (dx * dx + dy * dy < holeRadii[h] * holeRadii[h])
                {
                    line.x.push_back(std::numeric_limits<double>::quiet_NaN());
                    line.y.push_back(std::numeric_limits<double>::quiet_NaN());
                    return;
                }
            }
            line.x.push_back(x);
            line.y.push_back(y);
        };

        // Radial lines from origin to unit circle
        for (int i = 0; i < numRadials; ++i)
        {
            double angle = 2.0 * M_PI * i / numRadials;
            double cosA = std::cos(angle);
            double sinA = std::sin(angle);
            GridLine line;
            for (int j = 0; j <= pointsPerLine; ++j)
            {
                double r = static_cast<double>(j) / pointsPerLine;
                appendPoint(line, r * cosA, r * sinA);
            }
            m_sourceGridLines.push_back(std::move(line));
        }

        // Concentric circles at evenly spaced radii
        for (int i = 1; i <= numCircles; ++i)
        {
            double radius = static_cast<double>(i) / (numCircles + 1);
            GridLine line;
            for (int j = 0; j <= pointsPerLine; ++j)
            {
                double angle = 2.0 * M_PI * j / pointsPerLine;
                appendPoint(line, radius * std::cos(angle), radius * std::sin(angle));
            }
            m_sourceGridLines.push_back(std::move(line));
        }
        return;
    }

    // For domains not handled above, fall back to Cartesian grid with containment testing.
    BoundingBox bounds = computeBounds(m_sourceBoundaries);
    double marginX = (bounds.xMax - bounds.xMin) * 0.05;
    double marginY = (bounds.yMax - bounds.yMin) * 0.05;
    double xMin = bounds.xMin - marginX;
    double xMax = bounds.xMax + marginX;
    double yMin = bounds.yMin - marginY;
    double yMax = bounds.yMax + marginY;

    int pointsPerLine = 100;

    // Horizontal lines
    for (int i = 0; i < m_gridDensity; ++i)
    {
        bool containsFailureLogged = false;
        double y = yMin + (yMax - yMin) * (i + 0.5) / m_gridDensity;
        GridLine line;
        for (int j = 0; j <= pointsPerLine; ++j)
        {
            double x = xMin + (xMax - xMin) * static_cast<double>(j) / pointsPerLine;
            Complex z(x, y);
            bool inDomain = false;
            try
            {
                inDomain = sourceDomain->contains(z);
            }
            catch (const std::exception& e)
            {
                if (!containsFailureLogged)
                {
                    spdlog::warn("Domain contains() check failed at z=({}, {}): {}", x, y, e.what());
                    containsFailureLogged = true;
                }
            }

            if (inDomain)
            {
                line.x.push_back(x);
                line.y.push_back(y);
            }
            else
            {
                line.x.push_back(std::numeric_limits<double>::quiet_NaN());
                line.y.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        m_sourceGridLines.push_back(std::move(line));
    }

    // Vertical lines
    for (int i = 0; i < m_gridDensity; ++i)
    {
        bool containsFailureLogged = false;
        double x = xMin + (xMax - xMin) * (i + 0.5) / m_gridDensity;
        GridLine line;
        for (int j = 0; j <= pointsPerLine; ++j)
        {
            double y = yMin + (yMax - yMin) * static_cast<double>(j) / pointsPerLine;
            Complex z(x, y);
            bool inDomain = false;
            try
            {
                inDomain = sourceDomain->contains(z);
            }
            catch (const std::exception& e)
            {
                if (!containsFailureLogged)
                {
                    spdlog::warn("Domain contains() check failed at z=({}, {}): {}", x, y, e.what());
                    containsFailureLogged = true;
                }
            }

            if (inDomain)
            {
                line.x.push_back(x);
                line.y.push_back(y);
            }
            else
            {
                line.x.push_back(std::numeric_limits<double>::quiet_NaN());
                line.y.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        m_sourceGridLines.push_back(std::move(line));
    }
}

void VisualizationPanel::generateTargetGrid()
{
    m_targetGridLines.clear();

    if (!mp_currentMap)
    {
        return;
    }

    int totalPoints = 0;
    int failedPoints = 0;

    for (const auto& sourceLine : m_sourceGridLines)
    {
        GridLine targetLine;
        bool firstFailureLogged = false;
        for (size_t i = 0; i < sourceLine.x.size(); ++i)
        {
            double sx = sourceLine.x[i];
            double sy = sourceLine.y[i];

            if (std::isnan(sx) || std::isnan(sy))
            {
                targetLine.x.push_back(std::numeric_limits<double>::quiet_NaN());
                targetLine.y.push_back(std::numeric_limits<double>::quiet_NaN());
                continue;
            }

            ++totalPoints;
            Complex z(sx, sy);

            try
            {
                Complex w = mp_currentMap->map(z);
                targetLine.x.push_back(w.real());
                targetLine.y.push_back(w.imag());
            }
            catch (const std::invalid_argument& e)
            {
                ++failedPoints;
                targetLine.x.push_back(std::numeric_limits<double>::quiet_NaN());
                targetLine.y.push_back(std::numeric_limits<double>::quiet_NaN());
                if (!firstFailureLogged)
                {
                    spdlog::warn("Grid map evaluation configuration error: {}", e.what());
                    firstFailureLogged = true;
                }
            }
            catch (const std::runtime_error& e)
            {
                ++failedPoints;
                targetLine.x.push_back(std::numeric_limits<double>::quiet_NaN());
                targetLine.y.push_back(std::numeric_limits<double>::quiet_NaN());
                if (!firstFailureLogged)
                {
                    spdlog::debug("Grid map evaluation failed at z=({}, {}): {}", z.real(), z.imag(), e.what());
                    firstFailureLogged = true;
                }
            }
        }
        m_targetGridLines.push_back(std::move(targetLine));
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
            spdlog::debug("Grid generation: {}/{} points failed map evaluation", failedPoints, totalPoints);
        }
    }
}

// --- Rendering ---

void VisualizationPanel::renderSourceDomain()
{
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;

    std::string title = "Source Domain";
    if (mp_currentMap && mp_currentMap->getSourceDomain())
    {
        title = determinePlotTitle(mp_currentMap->getSourceDomain(), true);
    }

    if (ImPlot::BeginPlot(title.c_str(), plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");

        BoundingBox bounds = computeBounds(m_sourceBoundaries);
        double rangeX = bounds.xMax - bounds.xMin;
        double rangeY = bounds.yMax - bounds.yMin;
        double margin = std::max(rangeX, rangeY) * 0.15;
        ImPlot::SetupAxesLimits(bounds.xMin - margin, bounds.xMax + margin,
                                bounds.yMin - margin, bounds.yMax + margin,
                                ImGuiCond_FirstUseEver);

        // Plot grid lines first so boundaries render on top
        if (m_showGrid && !m_sourceGridLines.empty())
        {
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            plotGridLines(m_sourceGridLines);
            ImPlot::PopStyleColor();
        }

        // Plot boundary curves (on top of grid)
        for (const auto& curve : m_sourceBoundaries)
        {
            if (!curve.x.empty())
            {
                ImPlot::PlotLine(curve.label.c_str(),
                                 curve.x.data(), curve.y.data(),
                                 static_cast<int>(curve.x.size()));
            }
        }

        ImPlot::EndPlot();
    }
}

void VisualizationPanel::renderTargetDomain()
{
    ImVec2 plotSize = ImVec2(-1, -1);
    plotSize.x = ImGui::GetContentRegionAvail().x;

    std::string title = "Target Domain";
    if (mp_currentMap && mp_currentMap->getTargetDomain())
    {
        title = determinePlotTitle(mp_currentMap->getTargetDomain(), false);
    }

    if (ImPlot::BeginPlot(title.c_str(), plotSize, ImPlotFlags_Equal))
    {
        ImPlot::SetupAxes("Real", "Imaginary");

        BoundingBox bounds = computeBounds(m_targetBoundaries);
        double rangeX = bounds.xMax - bounds.xMin;
        double rangeY = bounds.yMax - bounds.yMin;
        double margin = std::max(rangeX, rangeY) * 0.15;
        ImPlot::SetupAxesLimits(bounds.xMin - margin, bounds.xMax + margin,
                                bounds.yMin - margin, bounds.yMax + margin,
                                ImGuiCond_FirstUseEver);

        // Plot mapped grid lines first so boundaries render on top
        if (m_showGrid && mp_currentMap && !m_targetGridLines.empty())
        {
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.7f, 0.7f, 0.7f, 0.8f));
            plotGridLines(m_targetGridLines);
            ImPlot::PopStyleColor();
        }

        // Plot boundary curves (on top of grid)
        for (const auto& curve : m_targetBoundaries)
        {
            if (!curve.x.empty())
            {
                ImPlot::PlotLine(curve.label.c_str(),
                                 curve.x.data(), curve.y.data(),
                                 static_cast<int>(curve.x.size()));
            }
        }

        // Show message if no map is loaded
        if (!mp_currentMap)
        {
            ImPlot::PlotText("No mapping computed", 0.0, 0.0);
        }

        ImPlot::EndPlot();
    }
}

void VisualizationPanel::clearGridData()
{
    m_sourceGridLines.clear();
    m_targetGridLines.clear();
}

void VisualizationPanel::plotGridLines(const std::vector<GridLine>& gridLines)
{
    // ImPlot's SkipNaN connected valid points across NaN gaps with a straight line,
    // which drew grid lines through holes and outside domain boundaries (GH-110).
    // Instead, split each grid line at NaN boundaries and plot each contiguous
    // segment as a separate PlotLine call.
    for (const auto& line : gridLines)
    {
        size_t n = line.x.size();
        size_t segStart = 0;

        while (segStart < n)
        {
            // Skip non-finite points (NaN from domain clipping, Inf from singular maps)
            while (segStart < n && !(std::isfinite(line.x[segStart]) && std::isfinite(line.y[segStart])))
            {
                ++segStart;
            }

            if (segStart >= n)
            {
                break;
            }

            // Find end of contiguous valid segment
            size_t segEnd = segStart;
            while (segEnd < n && std::isfinite(line.x[segEnd]) && std::isfinite(line.y[segEnd]))
            {
                ++segEnd;
            }

            // Need at least 2 points to draw a line segment
            int count = static_cast<int>(segEnd - segStart);
            if (count >= 2)
            {
                ImPlot::PlotLine("##grid", line.x.data() + segStart, line.y.data() + segStart, count);
            }

            segStart = segEnd;
        }
    }
}
