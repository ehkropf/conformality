#include "GuiController.h"
#include "VisualizationPanel.h"
#include "../core/ConformalMap.h"
#include "../domains/Domain.h"
#include "../methods/TheodorsenMethod.h"
#include <cmath>
#include <stdexcept>

GuiController::GuiController()
    : mp_currentMap{nullptr}
    , mp_sourceDomain{nullptr}
    , mp_targetDomain{nullptr}
    , mp_method{nullptr}
    , mp_visualizationPanel{nullptr}
    , m_ellipseA{2.0}
    , m_ellipseB{1.0}
    , m_samplePoints{256}
    , m_mappingType{MappingType::INTERIOR_TO_INTERIOR}
    , m_isComputing{false}
    , m_needsRecomputation{true}
    , m_lastComputationSuccessful{false}
    , m_lastConvergenceError{0.0}
    , m_lastErrorMessage{""}
{
}

GuiController::~GuiController()
{
    shutdown();
}

bool GuiController::initialize()
{
    // Create default domains and method
    createDomains();
    createMethod();
    
    return true;
}

void GuiController::shutdown()
{
    mp_currentMap.reset();
    mp_sourceDomain.reset();
    mp_targetDomain.reset();
    mp_method.reset();
    mp_visualizationPanel = nullptr;
}

void GuiController::setVisualizationPanel(VisualizationPanel* panel)
{
    mp_visualizationPanel = panel;
}

void GuiController::setEllipseParameters(double a, double b)
{
    if (a > 0.0 && b > 0.0 && (a != m_ellipseA || b != m_ellipseB))
    {
        m_ellipseA = a;
        m_ellipseB = b;
        m_needsRecomputation = true;
        
        // Update target domain
        createDomains();
        
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Ellipse parameters updated");
        }
    }
}

void GuiController::setSamplePoints(int points)
{
    if (isPowerOfTwo(points) && points != m_samplePoints)
    {
        m_samplePoints = points;
        m_needsRecomputation = true;
        
        // Recreate method with new sample points
        createMethod();
        
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Sample points updated to " + std::to_string(points));
        }
    }
}

void GuiController::setMappingType(MappingType type)
{
    if (type != m_mappingType)
    {
        m_mappingType = type;
        m_needsRecomputation = true;
        
        // Update domains based on mapping type
        createDomains();
        
        if (m_onStatusUpdate)
        {
            std::string typeStr = (type == MappingType::INTERIOR_TO_INTERIOR) ? "Internal" : "External";
            m_onStatusUpdate("Mapping type set to " + typeStr);
        }
    }
}

bool GuiController::computeMapping()
{
    if (m_isComputing)
    {
        return false;
    }
    
    if (!validateParameters())
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = "Invalid parameters";
        return false;
    }
    
    m_isComputing = true;
    m_lastComputationSuccessful = false;
    m_lastErrorMessage = "";
    
    if (m_onStatusUpdate)
    {
        m_onStatusUpdate("Computing conformal map...");
    }
    
    try
    {
        // Create fresh map with current parameters
        createMap();
        
        if (!mp_currentMap)
        {
            throw std::runtime_error("Failed to create conformal map");
        }
        
        // Compute the mapping (returns void, so we assume success if no exception)
        mp_currentMap->compute();
        bool success = true;
        
        if (success)
        {
            m_lastComputationSuccessful = true;
            m_needsRecomputation = false;
            
            // Update visualization
            updateVisualization();
            
            // Get convergence information from method
            if (auto theodorsenMethod = std::dynamic_pointer_cast<TheodorsenMethod>(mp_method))
            {
                // Note: This assumes TheodorsenMethod has a method to get convergence info
                // We'll need to add this to the TheodorsenMethod class
                m_lastConvergenceError = 0.001; // Placeholder
            }
            
            if (m_onStatusUpdate)
            {
                m_onStatusUpdate("Computation completed successfully");
            }
        }
        else
        {
            throw std::runtime_error("Computation failed to converge");
        }
    }
    catch (const std::exception& e)
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = e.what();
        
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Computation failed: " + m_lastErrorMessage);
        }
    }
    
    m_isComputing = false;
    
    if (m_onComputationComplete)
    {
        m_onComputationComplete();
    }
    
    return m_lastComputationSuccessful;
}

double GuiController::getEccentricity() const
{
    if (m_ellipseA >= m_ellipseB)
    {
        return std::sqrt(1.0 - (m_ellipseB * m_ellipseB) / (m_ellipseA * m_ellipseA));
    }
    else
    {
        return std::sqrt(1.0 - (m_ellipseA * m_ellipseA) / (m_ellipseB * m_ellipseB));
    }
}

void GuiController::setOnComputationComplete(std::function<void()> callback)
{
    m_onComputationComplete = callback;
}

void GuiController::setOnStatusUpdate(std::function<void(const std::string&)> callback)
{
    m_onStatusUpdate = callback;
}

void GuiController::createDomains()
{
    // Source domain is unit circle for Theodorsen method (what we're mapping FROM)
    mp_sourceDomain = std::make_shared<CircularDomain>(Complex(0.0, 0.0), 1.0, false);
    
    // Target domain is ellipse (what we're mapping TO)
    mp_targetDomain = std::make_shared<EllipticalDomain>(
        m_ellipseA,
        m_ellipseB,
        0.0,  // No rotation for now
        Complex(0.0, 0.0),
        m_mappingType == MappingType::EXTERIOR_TO_INTERIOR  // External domain for external mapping
    );
}

void GuiController::createMethod()
{
    mp_method = std::make_shared<TheodorsenMethod>(m_samplePoints);
}

void GuiController::createMap()
{
    if (!mp_sourceDomain || !mp_targetDomain || !mp_method)
    {
        createDomains();
        createMethod();
    }
    
    mp_currentMap = std::make_shared<ConformalMap>(
        mp_sourceDomain,
        mp_targetDomain,
        mp_method
    );
}

void GuiController::updateVisualization()
{
    if (mp_visualizationPanel && mp_currentMap)
    {
        mp_visualizationPanel->updateMap(mp_currentMap);
    }
}

bool GuiController::validateParameters()
{
    // Check ellipse parameters
    if (m_ellipseA <= 0.0 || m_ellipseB <= 0.0)
    {
        return false;
    }
    
    // Check sample points is power of 2
    if (!isPowerOfTwo(m_samplePoints))
    {
        return false;
    }
    
    // Check minimum sample points
    if (m_samplePoints < 16)
    {
        return false;
    }
    
    return true;
}

bool GuiController::isPowerOfTwo(int n) const
{
    return n > 0 && (n & (n - 1)) == 0;
}