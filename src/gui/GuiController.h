#pragma once

#include <memory>
#include <functional>
#include <string>

class ConformalMap;
class VisualizationPanel;

/**
 * @brief Controller class that bridges GUI components with computational backend
 *
 * Manages the lifecycle of conformal maps, coordinates between GUI controls
 * and computation, and handles parameter updates and map generation.
 *
 * The controller starts in a blank state and accepts any ConformalMap via
 * the loadMap() method. It is method-agnostic.
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

    // Computation state
    bool m_isComputing;

    // Computation results
    bool m_lastComputationSuccessful;
    double m_lastConvergenceError;
    std::string m_lastErrorMessage;

    // Callbacks for GUI updates
    std::function<void()> m_onComputationComplete;
    std::function<void(const std::string&)> m_onStatusUpdate;

public:
    GuiController();
    ~GuiController();

    bool initialize();
    void shutdown();

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
     * @brief Clear the current map and reset state
     */
    void clear();

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
     * @brief Trigger computation of conformal map
     * @return true if computation started successfully
     */
    bool computeMapping();

    /**
     * @brief Check if computation is currently running
     * @return true if computing
     */
    bool isComputing() const { return m_isComputing; }

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
     * @brief Set callback for computation completion
     * @param callback Function to call when computation completes
     */
    void setOnComputationComplete(std::function<void()> callback);

    /**
     * @brief Set callback for status updates
     * @param callback Function to call with status messages
     */
    void setOnStatusUpdate(std::function<void(const std::string&)> callback);

private:
    void updateVisualization();
};
