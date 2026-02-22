#include "GuiController.h"
#include "VisualizationPanel.h"
#include "../core/ConformalMap.h"

GuiController::GuiController()
    : mp_currentMap{nullptr}
    , mp_visualizationPanel{nullptr}
    , m_isComputing{false}
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
    return true;
}

void GuiController::shutdown()
{
    mp_currentMap.reset();
    m_mapDescription.clear();
    mp_visualizationPanel = nullptr;
}

void GuiController::setVisualizationPanel(VisualizationPanel* panel)
{
    mp_visualizationPanel = panel;
}

void GuiController::loadMap(std::shared_ptr<ConformalMap> map, const std::string& description)
{
    mp_currentMap = map;
    m_mapDescription = description;
    m_lastComputationSuccessful = false;
    m_lastConvergenceError = 0.0;
    m_lastErrorMessage.clear();

    updateVisualization();

    if (m_onStatusUpdate)
    {
        m_onStatusUpdate(description.empty() ? "Map loaded" : description);
    }
}

void GuiController::clear()
{
    mp_currentMap.reset();
    m_mapDescription.clear();
    m_lastComputationSuccessful = false;
    m_lastConvergenceError = 0.0;
    m_lastErrorMessage.clear();

    if (m_onStatusUpdate)
    {
        m_onStatusUpdate("Ready");
    }
}

bool GuiController::hasMap() const
{
    return mp_currentMap != nullptr;
}

std::shared_ptr<ConformalMap> GuiController::getCurrentMap() const
{
    return mp_currentMap;
}

const std::string& GuiController::getMapDescription() const
{
    return m_mapDescription;
}

bool GuiController::computeMapping()
{
    if (!mp_currentMap)
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = "No map loaded";
        return false;
    }

    if (m_isComputing)
    {
        return false;
    }

    m_isComputing = true;
    m_lastComputationSuccessful = false;
    m_lastConvergenceError = 0.0;
    m_lastErrorMessage.clear();

    if (m_onStatusUpdate)
    {
        m_onStatusUpdate("Computing conformal map...");
    }

    try
    {
        mp_currentMap->compute();
        m_lastComputationSuccessful = true;

        updateVisualization();

        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Computation completed successfully");
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

void GuiController::setOnComputationComplete(std::function<void()> callback)
{
    m_onComputationComplete = callback;
}

void GuiController::setOnStatusUpdate(std::function<void(const std::string&)> callback)
{
    m_onStatusUpdate = callback;
}

void GuiController::updateVisualization()
{
    if (mp_visualizationPanel && mp_currentMap)
    {
        mp_visualizationPanel->updateMap(mp_currentMap);
    }
}
