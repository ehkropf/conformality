#pragma once

#include <memory>
#include <string>

class VisualizationPanel;
class GuiController;
class Application;

class MainWindow
{
private:
    bool m_showDemoWindow;
    bool m_showAboutDialog;
    std::unique_ptr<VisualizationPanel> mp_visualizationPanel;
    std::unique_ptr<GuiController> mp_controller;
    std::string m_statusMessage;
    std::string m_computationPhase;
    Application* mp_application;

public:
    MainWindow();
    ~MainWindow();

    bool initialize();
    void render();
    void shutdown();

    void setApplication(Application* app) { mp_application = app; }

private:
    void renderMenuBar();
    void renderMainLayout();
    void renderControlPanel();
    void renderVisualizationPanel();
    void renderStatusPanel();
    void renderAboutDialog();

    // GUI state for controls
    bool m_showGrid;
    int m_gridDensity;

    // Cached per-frame live progress (read once in renderMainLayout)
    int m_cachedLiveIter{0};
    double m_cachedLiveResidual{0.0};

    // Callback for status updates
    void onStatusUpdate(const std::string& message);
    void onComputationComplete();
};
