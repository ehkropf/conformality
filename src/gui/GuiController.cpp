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
    cancelAndJoin();
    mp_currentMap.reset();
    m_mapDescription.clear();
    mp_visualizationPanel = nullptr;
}

void GuiController::update()
{
    // Drain message queue → fire status callbacks on GUI thread
    pollStatusMessages();

    // Check if background computation finished
    if (!m_isComputing.load() && m_computeThread.joinable())
    {
        m_computeThread.join();
        updateVisualization();

        if (m_onComputationComplete)
        {
            m_onComputationComplete();
        }
    }
}

void GuiController::setVisualizationPanel(VisualizationPanel* panel)
{
    mp_visualizationPanel = panel;
}

void GuiController::loadMap(std::shared_ptr<ConformalMap> map, const std::string& description)
{
    cancelAndJoin();

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
    cancelAndJoin();

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

void GuiController::cancelComputation()
{
    m_cancelRequested.store(true);
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

    if (m_isComputing.load())
    {
        return false;
    }

    // Reset result state before launching
    m_lastComputationSuccessful = false;
    m_lastConvergenceError = 0.0;
    m_lastErrorMessage.clear();
    m_lastIterationCount = 0;
    m_hasConverged = false;
    m_cancelRequested.store(false);
    {
        std::lock_guard<std::mutex> lock(m_progressMutex);
        m_liveProgress = ComputationProgress{};
    }

    // Join any previous thread that hasn't been joined yet
    if (m_computeThread.joinable())
    {
        m_computeThread.join();
    }

    m_isComputing.store(true);
    postStatusMessage("Computing conformal map...");

    m_computeThread = std::thread(&GuiController::computeInBackground, this);

    return true;
}

void GuiController::getLiveProgress(int& iteration, double& residual) const
{
    std::lock_guard<std::mutex> lock(m_progressMutex);
    iteration = m_liveProgress.currentIteration;
    residual = m_liveProgress.currentResidual;
}

void GuiController::computeInBackground()
{
    // Reset live progress
    {
        std::lock_guard<std::mutex> lock(m_progressMutex);
        m_liveProgress = ComputationProgress{};
    }

    try
    {
        // Wire cancellation check into FornbergMC if applicable
        auto fm = std::dynamic_pointer_cast<FornbergMC>(mp_currentMap->getMethod());
        if (fm)
        {
            fm->setCancellationCheck([this]() { return m_cancelRequested.load(); });
        }

        // Wire StatusManager callback for live progress updates
        auto status_manager = fm ? std::dynamic_pointer_cast<StatusManager>(fm->getStatusManager()) : nullptr;
        if (status_manager)
        {
            status_manager->setStatusCallback([this](const StatusMessage& msg) {
                postStatusMessage("[" + msg.component + "] " + msg.message);

                // Parse iteration progress from FornbergMC INFO messages
                if (msg.component == "FornbergMC" && msg.level == StatusLevel::INFO
                    && msg.message.find("Newton iteration") != std::string::npos)
                {
                    std::lock_guard<std::mutex> lock(m_progressMutex);
                    // Parse "Newton iteration N: residual=X"
                    auto colon_pos = msg.message.find(':');
                    auto eq_pos = msg.message.find("residual=");
                    if (colon_pos != std::string::npos)
                    {
                        try
                        {
                            auto iter_str = msg.message.substr(17, colon_pos - 17);
                            m_liveProgress.currentIteration = std::stoi(iter_str);
                        }
                        catch (...) {}
                    }
                    if (eq_pos != std::string::npos)
                    {
                        try
                        {
                            m_liveProgress.currentResidual = std::stod(msg.message.substr(eq_pos + 9));
                        }
                        catch (...) {}
                    }
                }
            });
        }

        mp_currentMap->compute();

        // Extract results (still on worker thread, but before m_isComputing goes false)
        m_lastComputationSuccessful = true;
        if (fm)
        {
            m_lastConvergenceError = fm->getCurrentResidual();
            m_lastIterationCount = static_cast<int>(fm->getResidualHistory().size());
            m_hasConverged = fm->hasConverged();
        }

        postStatusMessage("Computation completed successfully");
    }
    catch (const std::invalid_argument& e)
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = e.what();
        postStatusMessage("Configuration error: " + m_lastErrorMessage);
    }
    catch (const std::runtime_error& e)
    {
        m_lastComputationSuccessful = false;
        m_lastErrorMessage = e.what();
        postStatusMessage("Computation failed: " + m_lastErrorMessage);
    }

    // Clear cancellation check and StatusManager callback to avoid dangling references
    auto fm_cleanup = std::dynamic_pointer_cast<FornbergMC>(mp_currentMap->getMethod());
    if (fm_cleanup)
    {
        fm_cleanup->setCancellationCheck(nullptr);
        auto sm = std::dynamic_pointer_cast<StatusManager>(fm_cleanup->getStatusManager());
        if (sm)
        {
            sm->setStatusCallback(nullptr);
        }
    }

    // Release happens-before: GUI thread reads results after observing m_isComputing == false
    m_isComputing.store(false);
}

void GuiController::postStatusMessage(const std::string& message)
{
    std::lock_guard<std::mutex> lock(m_messageQueueMutex);
    m_messageQueue.push_back(message);
}

void GuiController::pollStatusMessages()
{
    std::deque<std::string> messages;
    {
        std::lock_guard<std::mutex> lock(m_messageQueueMutex);
        messages.swap(m_messageQueue);
    }

    for (const auto& msg : messages)
    {
        if (m_onStatusUpdate)
        {
            m_onStatusUpdate(msg);
        }
    }
}

void GuiController::cancelAndJoin()
{
    if (m_isComputing.load())
    {
        m_cancelRequested.store(true);
    }
    if (m_computeThread.joinable())
    {
        m_computeThread.join();
    }
    m_cancelRequested.store(false);
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
