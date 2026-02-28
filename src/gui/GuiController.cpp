#include "GuiController.h"
#include "VisualizationPanel.h"
#include "../core/ConformalMap.h"
#include "../core/StatusManager.h"
#include "../domains/FornbergCanonicalDomain.h"
#include "../examples/ThesisExamples.h"
#include "../methods/FornbergMC.h"
#include "../methods/PMatrixBuilder.h"
#include "../numerics/CGSolver.h"

GuiController::GuiController()
    : mp_currentMap{nullptr}
    , mp_visualizationPanel{nullptr}
    , m_isComputing{false}
    , m_lastComputationSuccessful{false}
    , m_lastConvergenceError{0.0}
    , m_lastErrorMessage{""}
    , m_lastIterationCount{0}
    , m_hasConverged{false}
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
    m_lastIterationCount = 0;
    m_hasConverged = false;

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
    m_lastIterationCount = 0;
    m_hasConverged = false;

    if (m_onStatusUpdate)
    {
        m_onStatusUpdate("Ready");
    }
}

void GuiController::loadThesisExample(int exampleNumber)
{
    using namespace conformality::examples;

    clear();

    try
    {
        auto preset = ThesisExamples::getExample(exampleNumber);

        auto source_domain = std::make_shared<FornbergCanonicalDomain>(
            preset.initial_centers, preset.initial_radii, preset.config.N);

        auto method = std::make_shared<FornbergMC>(preset.config);
        auto status_manager = std::make_shared<StatusManager>();
        status_manager->enableLogging(LogOutput::CONSOLE);
        method->setStatusManager(status_manager);

        auto map = std::make_shared<ConformalMap>(source_domain, preset.target_domain, method);
        loadMap(map, preset.name + ": " + preset.description);
    }
    catch (const std::invalid_argument& e)
    {
        m_lastErrorMessage = e.what();
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Failed to load example: " + m_lastErrorMessage);
        }
    }
    catch (const std::runtime_error& e)
    {
        m_lastErrorMessage = e.what();
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Failed to load example: " + m_lastErrorMessage);
        }
    }
}

void GuiController::reset()
{
    clear();

    if (mp_visualizationPanel)
    {
        mp_visualizationPanel->updateMap(nullptr);
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
    m_lastIterationCount = 0;
    m_hasConverged = false;

    if (m_onStatusUpdate)
    {
        m_onStatusUpdate("Computing conformal map...");
    }

    try
    {
        mp_currentMap->compute();
        m_lastComputationSuccessful = true;

        auto fm = std::dynamic_pointer_cast<FornbergMC>(mp_currentMap->getMethod());
        if (fm)
        {
            m_lastConvergenceError = fm->getCurrentResidual();
            m_lastIterationCount = static_cast<int>(fm->getResidualHistory().size());
            m_hasConverged = fm->hasConverged();
        }

        updateVisualization();

        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Computation completed successfully");
        }
    }
    catch (const std::invalid_argument& e)
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = e.what();

        if (m_onStatusUpdate)
        {
            m_onStatusUpdate("Configuration error: " + m_lastErrorMessage);
        }
    }
    catch (const std::runtime_error& e)
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
