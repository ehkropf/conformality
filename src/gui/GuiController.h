#pragma once

#include "../core/MethodInfo.h"

#include <atomic>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>

class ConformalMap;
class VisualizationPanel;

/**
 * @brief Controller class that bridges GUI components with computational backend
 *
 * Manages the lifecycle of conformal maps, coordinates between GUI controls
 * and computation, and handles parameter updates and map generation.
 *
 * The controller starts in a blank state and accepts any ConformalMap via
 * the loadMap() method. After computation, it extracts results via the
 * MethodInfo interface when available.
 *
 * Computation runs on a background thread to keep the GUI responsive.
 * Call update() each frame to process status messages and detect completion.
 */
class GuiController
{
private:
    // Current conformal mapping
    std::shared_ptr<ConformalMap> mp_currentMap;

    // GUI components
    VisualizationPanel* mp_visualizationPanel;

    // Map description
    std::string m_mapDescription;

    // Computation state (atomic for cross-thread access)
    std::atomic<bool> m_isComputing{false};
    std::atomic<bool> m_cancelRequested{false};
    std::thread m_computeThread;

    // Thread-safe message queue for status updates from worker thread
    std::mutex m_messageQueueMutex;
    std::deque<std::string> m_messageQueue;

    // Live progress from worker thread (protected by m_progressMutex)
    struct ComputationProgress
    {
        int currentIteration{0};    ///< 1-based Newton iteration number (0 = not started)
        double currentResidual{0.0}; ///< Newton update residual norm
    };
    mutable std::mutex m_progressMutex;
    ComputationProgress m_liveProgress;

    // Computation results (written by worker thread before m_isComputing goes false)
    bool m_lastComputationSuccessful{false};
    double m_lastConvergenceError{0.0};
    std::string m_lastErrorMessage;
    int m_lastIterationCount{0};
    bool m_hasConverged{false};
    conformality::MethodInfo m_lastMethodInfo;

    // Callbacks for GUI updates (only called from main/GUI thread)
    std::function<void()> m_onComputationComplete;
    std::function<void(const std::string&)> m_onStatusUpdate;

    void postStatusMessage(const std::string& message);
    void pollStatusMessages();
    void computeInBackground();
    void cancelAndJoin();
    void updateVisualization();

public:
    GuiController();
    ~GuiController();

    bool initialize();
    void shutdown();

    /**
     * @brief Poll for background thread completion and drain message queue
     *
     * Must be called each frame from the render loop. Delivers queued status
     * messages to the onStatusUpdate callback, and when computation finishes,
     * joins the thread, updates visualization, and fires onComputationComplete.
     */
    void update();

    /**
     * @brief Set the visualization panel to update
     * @param panel Pointer to visualization panel
     */
    void setVisualizationPanel(VisualizationPanel* panel);

    /**
     * @brief Load a conformal map into the controller
     * @param map The conformal map to load
     * @param description Optional human-readable description
     */
    void loadMap(std::shared_ptr<ConformalMap> map, const std::string& description = "");

    /**
     * @brief Load a thesis example by number and prepare it for computation
     * @param exampleNumber Thesis example number (see ThesisExamples::availableExamples())
     */
    void loadThesisExample(int exampleNumber);

    /**
     * @brief Full user-facing reset: clears map, state, and visualization
     */
    void reset();

    /**
     * @brief Clear map and computation state (use reset() to also clear visualization)
     */
    void clear();

    /**
     * @brief Request cancellation of a running computation
     *
     * Sets a flag that the worker thread checks cooperatively. The computation
     * will stop at the next cancellation check point (top of Newton loop).
     */
    void cancelComputation();

    /**
     * @brief Check if a map is currently loaded
     * @return true if a map is loaded
     */
    bool hasMap() const;

    /**
     * @brief Get the currently loaded map
     * @return Shared pointer to current map (may be null)
     */
    std::shared_ptr<ConformalMap> getCurrentMap() const;

    /**
     * @brief Get the description of the currently loaded map
     * @return Map description string
     */
    const std::string& getMapDescription() const;

    /**
     * @brief Launch background computation of conformal map
     * @return true if computation was started, false if already computing or no map
     */
    bool computeMapping();

    /**
     * @brief Check if computation is currently running
     * @return true if computing
     */
    bool isComputing() const { return m_isComputing.load(); }

    /**
     * @brief Get a snapshot of live computation progress (thread-safe)
     * @param iteration Output: 1-based Newton iteration number (0 if no iteration reported yet)
     * @param residual Output: current Newton residual (meaningful only when iteration > 0)
     */
    void getLiveProgress(int& iteration, double& residual) const;

    /**
     * @brief Check if last computation was successful
     * @return true if successful
     */
    bool wasLastComputationSuccessful() const { return m_lastComputationSuccessful; }

    /**
     * @brief Get last convergence error
     * @return Convergence error value
     */
    double getLastConvergenceError() const { return m_lastConvergenceError; }

    /**
     * @brief Get last error message
     * @return Error message string
     */
    const std::string& getLastErrorMessage() const { return m_lastErrorMessage; }

    /**
     * @brief Get iteration count from last computation
     * @return Number of iterations (0 if not yet computed)
     */
    int getLastIterationCount() const { return m_lastIterationCount; }

    /**
     * @brief Check if last computation converged
     * @return true if converged within tolerance
     */
    bool hasConverged() const { return m_hasConverged; }

    /**
     * @brief Check if last computation was cancelled by the user
     * @return true if cancelled
     */
    bool wasCancelled() const { return m_cancelRequested.load(); }

    /**
     * @brief Get method info from last computation
     * @return conformality::MethodInfo with parameters and results
     */
    const conformality::MethodInfo& getLastMethodInfo() const { return m_lastMethodInfo; }

    /**
     * @brief Set callback for computation completion
     * @param callback Function to call when computation completes
     */
    void setOnComputationComplete(std::function<void()> callback);

    /**
     * @brief Set callback for status updates
     * @param callback Function to call with status messages
     */
    void setOnStatusUpdate(std::function<void(const std::string&)> callback);
};
