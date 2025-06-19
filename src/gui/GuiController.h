#pragma once

#include "../core/Types.h"
#include <memory>
#include <functional>

class ConformalMap;
class Domain;
class ConformalMapMethod;
class VisualizationPanel;

/**
 * @brief Controller class that bridges GUI components with computational backend
 * 
 * Manages the lifecycle of conformal maps, coordinates between GUI controls
 * and computation, and handles parameter updates and map generation.
 */
class GuiController
{
private:
    // Current conformal mapping setup
    std::shared_ptr<ConformalMap> mp_currentMap;
    std::shared_ptr<Domain> mp_sourceDomain;
    std::shared_ptr<Domain> mp_targetDomain;
    std::shared_ptr<ConformalMapMethod> mp_method;
    
    // GUI components
    VisualizationPanel* mp_visualizationPanel;
    
    // Current parameters
    double m_ellipseA;
    double m_ellipseB;
    int m_samplePoints;
    MappingType m_mappingType;
    bool m_isComputing;
    bool m_needsRecomputation;
    
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
     * @brief Update ellipse parameters
     * @param a Semi-major axis
     * @param b Semi-minor axis
     */
    void setEllipseParameters(double a, double b);
    
    /**
     * @brief Set number of sample points (must be power of 2)
     * @param points Number of sample points
     */
    void setSamplePoints(int points);
    
    /**
     * @brief Set mapping type (internal or external)
     * @param type Mapping type
     */
    void setMappingType(MappingType type);
    
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
     * @brief Get current ellipse semi-major axis
     * @return Semi-major axis value
     */
    double getEllipseA() const { return m_ellipseA; }
    
    /**
     * @brief Get current ellipse semi-minor axis
     * @return Semi-minor axis value
     */
    double getEllipseB() const { return m_ellipseB; }
    
    /**
     * @brief Get current number of sample points
     * @return Number of sample points
     */
    int getSamplePoints() const { return m_samplePoints; }
    
    /**
     * @brief Get current mapping type
     * @return Mapping type
     */
    MappingType getMappingType() const { return m_mappingType; }
    
    /**
     * @brief Get eccentricity of current ellipse
     * @return Eccentricity value
     */
    double getEccentricity() const;
    
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
    void createDomains();
    void createMethod();
    void createMap();
    void updateVisualization();
    bool validateParameters();
    
    // Helper to check if number is power of 2
    bool isPowerOfTwo(int n) const;
};