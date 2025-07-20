#include "MainWindow.h"
#include "VisualizationPanel.h"
#include "GuiController.h"
#include "Application.h"
#include "../core/Types.h"
#include "imgui.h"
#include "implot.h"
#include <cmath>

MainWindow::MainWindow()
    : m_showDemoWindow{false}
    , m_showAboutDialog{false}
    , mp_visualizationPanel{nullptr}
    , mp_controller{nullptr}
    , m_statusMessage{"Ready"}
    , mp_application{nullptr}
    , m_ellipseA{2.0f}
    , m_ellipseB{1.0f}
    , m_samplePoints{256}
    , m_mappingTypeIndex{0}
    , m_showGrid{true}
{
}

MainWindow::~MainWindow()
{
    shutdown();
}

bool MainWindow::initialize()
{
    // Initialize ImPlot
    ImPlot::CreateContext();
    
    // Create visualization panel
    mp_visualizationPanel = std::make_unique<VisualizationPanel>();
    if (!mp_visualizationPanel->initialize())
    {
        return false;
    }
    
    // Create controller
    mp_controller = std::make_unique<GuiController>();
    if (!mp_controller->initialize())
    {
        return false;
    }
    
    // Connect visualization panel to controller
    mp_controller->setVisualizationPanel(mp_visualizationPanel.get());
    
    // Set up callbacks
    mp_controller->setOnStatusUpdate([this](const std::string& message) {
        onStatusUpdate(message);
    });
    
    mp_controller->setOnComputationComplete([this]() {
        onComputationComplete();
    });
    
    // Set initial parameters
    mp_controller->setEllipseParameters(m_ellipseA, m_ellipseB);
    mp_controller->setSamplePoints(m_samplePoints);
    mp_controller->setMappingType(m_mappingTypeIndex == 0 ? MappingType::INTERIOR_TO_INTERIOR : MappingType::EXTERIOR_TO_INTERIOR);
    
    return true;
}

void MainWindow::render()
{
    renderMenuBar();
    renderMainLayout();
    
    // Show ImGui demo window for development/testing
    if (m_showDemoWindow)
    {
        ImGui::ShowDemoWindow(&m_showDemoWindow);
    }
    
    // Show About dialog
    if (m_showAboutDialog)
    {
        renderAboutDialog();
    }
}

void MainWindow::shutdown()
{
    if (mp_controller)
    {
        mp_controller->shutdown();
        mp_controller.reset();
    }
    
    if (mp_visualizationPanel)
    {
        mp_visualizationPanel->shutdown();
        mp_visualizationPanel.reset();
    }
    
    // Cleanup ImPlot context
    ImPlot::DestroyContext();
}

void MainWindow::renderMenuBar()
{
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("New"))
            {
                // Future: Reset to default configuration
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit"))
            {
                if (mp_application)
                {
                    mp_application->requestClose();
                }
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("View"))
        {
            if (ImGui::MenuItem("Show Grid", nullptr, &m_showGrid))
            {
                if (mp_visualizationPanel)
                {
                    mp_visualizationPanel->setGridVisible(m_showGrid);
                }
            }
            ImGui::Separator();
            ImGui::MenuItem("Demo Window", nullptr, &m_showDemoWindow);
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("Help"))
        {
            if (ImGui::MenuItem("About"))
            {
                m_showAboutDialog = true;
            }
            ImGui::EndMenu();
        }
        
        ImGui::EndMainMenuBar();
    }
}

void MainWindow::renderMainLayout()
{
    // Get the main viewport
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(viewport->Pos.x, viewport->Pos.y + ImGui::GetFrameHeight()));
    ImGui::SetNextWindowSize(ImVec2(viewport->Size.x, viewport->Size.y - ImGui::GetFrameHeight()));
    
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | 
                                   ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings;
    
    if (ImGui::Begin("MainLayout", nullptr, window_flags))
    {
        // Three-panel layout: Control | Visualization | Status
        
        // Control Panel (left column)
        renderControlPanel();
        
        ImGui::SameLine();
        
        // Visualization Panel (center)
        renderVisualizationPanel();
        
        // Status Panel (bottom)
        renderStatusPanel();
    }
    ImGui::End();
}

void MainWindow::renderControlPanel()
{
    float control_width = 300.0f;
    float status_height = 100.0f;
    
    if (ImGui::BeginChild("ControlPanel", ImVec2(control_width, -status_height), true))
    {
        ImGui::Text("Control Panel");
        ImGui::Separator();
        
        // Target ellipse parameters
        bool parametersChanged = false;
        
        if (ImGui::SliderFloat("Target Ellipse a", &m_ellipseA, 0.1f, 10.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp))
        {
            parametersChanged = true;
        }
        
        if (ImGui::SliderFloat("Target Ellipse b", &m_ellipseB, 0.1f, m_ellipseA, "%.2f", ImGuiSliderFlags_AlwaysClamp))
        {
            parametersChanged = true;
        }
        
        // Show eccentricity
        double eccentricity = 0.0;
        if (mp_controller)
        {
            eccentricity = mp_controller->getEccentricity();
        }
        ImGui::Text("Eccentricity: %.3f", eccentricity);
        
        ImGui::Separator();
        
        // Sample points (powers of 2)
        const char* sampleOptions[] = {"16", "32", "64", "128", "256", "512", "1024", "2048"};
        int sampleValues[] = {16, 32, 64, 128, 256, 512, 1024, 2048};
        int currentSampleIndex = 4; // Default to 256
        
        // Find current index
        for (int i = 0; i < 8; ++i)
        {
            if (sampleValues[i] == m_samplePoints)
            {
                currentSampleIndex = i;
                break;
            }
        }
        
        if (ImGui::Combo("Sample Points", &currentSampleIndex, sampleOptions, 8))
        {
            m_samplePoints = sampleValues[currentSampleIndex];
            if (mp_controller)
            {
                mp_controller->setSamplePoints(m_samplePoints);
            }
        }
        
        // Mapping type
        const char* mappingTypes[] = {"Internal", "External"};
        if (ImGui::Combo("Mapping Type", &m_mappingTypeIndex, mappingTypes, 2))
        {
            if (mp_controller)
            {
                MappingType type = (m_mappingTypeIndex == 0) ? MappingType::INTERIOR_TO_INTERIOR : MappingType::EXTERIOR_TO_INTERIOR;
                mp_controller->setMappingType(type);
            }
        }
        
        ImGui::Separator();
        
        // Grid controls
        if (ImGui::Checkbox("Show Grid", &m_showGrid))
        {
            if (mp_visualizationPanel)
            {
                mp_visualizationPanel->setGridVisible(m_showGrid);
            }
        }
        
        // Grid density
        static int gridDensity = 8;
        if (ImGui::SliderInt("Grid Density", &gridDensity, 4, 16))
        {
            if (mp_visualizationPanel)
            {
                mp_visualizationPanel->setGridDensity(gridDensity);
            }
        }
        
        ImGui::Separator();
        
        // Compute button
        bool isComputing = mp_controller ? mp_controller->isComputing() : false;
        
        if (isComputing)
        {
            ImGui::Text("Computing...");
        }
        else
        {
            if (ImGui::Button("Compute Mapping", ImVec2(-1, 0)))
            {
                if (mp_controller)
                {
                    mp_controller->computeMapping();
                }
            }
        }
        
        // Update parameters if they changed
        if (parametersChanged && mp_controller)
        {
            mp_controller->setEllipseParameters(m_ellipseA, m_ellipseB);
        }
    }
    ImGui::EndChild();
}

void MainWindow::renderVisualizationPanel()
{
    float status_height = 100.0f;
    
    if (ImGui::BeginChild("VisualizationPanel", ImVec2(0, -status_height), true))
    {
        if (mp_visualizationPanel)
        {
            mp_visualizationPanel->render();
        }
        else
        {
            ImGui::Text("Visualization Panel");
            ImGui::Separator();
            ImGui::Text("Loading...");
        }
    }
    ImGui::EndChild();
}

void MainWindow::renderStatusPanel()
{
    if (ImGui::BeginChild("StatusPanel", ImVec2(0, 0), true))
    {
        ImGui::Text("Status: %s", m_statusMessage.c_str());
        
        if (mp_controller)
        {
            ImGui::SameLine();
            if (mp_controller->wasLastComputationSuccessful())
            {
                ImGui::Text("| Convergence: %.6f", mp_controller->getLastConvergenceError());
            }
            else if (!mp_controller->getLastErrorMessage().empty())
            {
                ImGui::Text("| Error: %s", mp_controller->getLastErrorMessage().c_str());
            }
        }
        
        ImGui::SameLine();
        ImGui::Text("| Performance: %.3f ms/frame (%.1f FPS)", 
                   1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    }
    ImGui::EndChild();
}

void MainWindow::onStatusUpdate(const std::string& message)
{
    m_statusMessage = message;
}

void MainWindow::onComputationComplete()
{
    // Additional actions when computation completes can be added here
}

void MainWindow::renderAboutDialog()
{
    if (m_showAboutDialog)
    {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
        
        if (ImGui::Begin("About Conformality Mapping Tool", &m_showAboutDialog, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize))
        {
            ImGui::Text("Conformality Mapping Tool");
            ImGui::Text("Version: 1/∞");
            ImGui::Separator();
            ImGui::Text("A tool for visualizing conformal mappings");
            ImGui::Text("from the unit circle to elliptical domains.");
            ImGui::Separator();
            
            if (ImGui::Button("Close", ImVec2(120, 0)))
            {
                m_showAboutDialog = false;
            }
        }
        ImGui::End();
    }
}