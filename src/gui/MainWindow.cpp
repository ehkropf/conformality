#include "MainWindow.h"
#include "imgui.h"
#include <cmath>

MainWindow::MainWindow()
    : m_showDemoWindow{false}
{
}

MainWindow::~MainWindow()
{
    shutdown();
}

bool MainWindow::initialize()
{
    // Panel initialization will be added when implementing individual panels
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
}

void MainWindow::shutdown()
{
    // Since panels are forward declared and not yet implemented,
    // they will be nullptr. The unique_ptr destructor will handle cleanup
    // when the actual panel classes are implemented.
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
                // Future: Signal application to close
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("View"))
        {
            if (ImGui::MenuItem("Show Grid"))
            {
                // Future: Toggle grid visibility
            }
            ImGui::Separator();
            ImGui::MenuItem("Demo Window", nullptr, &m_showDemoWindow);
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("Help"))
        {
            if (ImGui::MenuItem("About"))
            {
                // Future: Show about dialog
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
        // Create a simple test layout for now
        ImGui::Text("Conformality Mapping Tool");
        ImGui::Separator();
        
        // Three-panel layout: Control | Visualization | Status
        float control_width = 300.0f;
        float status_height = 100.0f;
        
        // Control Panel (left column)
        if (ImGui::BeginChild("ControlPanel", ImVec2(control_width, -status_height), true))
        {
            ImGui::Text("Control Panel");
            ImGui::Separator();
            
            // Placeholder controls
            static float ellipse_a = 2.0f;
            static float ellipse_b = 1.0f;
            
            ImGui::SliderFloat("Ellipse a", &ellipse_a, 1.0f, 5.0f);
            ImGui::SliderFloat("Ellipse b", &ellipse_b, 0.1f, ellipse_a);
            
            float eccentricity = sqrt(1.0f - (ellipse_b * ellipse_b) / (ellipse_a * ellipse_a));
            ImGui::Text("Eccentricity: %.3f", eccentricity);
            
            ImGui::Separator();
            if (ImGui::Button("Compute Mapping"))
            {
                // Future: Trigger computation
            }
        }
        ImGui::EndChild();
        
        ImGui::SameLine();
        
        // Visualization Panel (center)
        if (ImGui::BeginChild("VisualizationPanel", ImVec2(0, -status_height), true))
        {
            ImGui::Text("Visualization Panel");
            ImGui::Separator();
            ImGui::Text("Canonical Domain | Target Domain");
            ImGui::Text("(Plots will go here)");
        }
        ImGui::EndChild();
        
        // Status Panel (bottom)
        if (ImGui::BeginChild("StatusPanel", ImVec2(0, 0), true))
        {
            ImGui::Text("Status: Ready");
            ImGui::SameLine();
            ImGui::Text("| Convergence: N/A");
            ImGui::SameLine();
            ImGui::Text("| Performance: %.3f ms/frame (%.1f FPS)", 
                       1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        }
        ImGui::EndChild();
    }
    ImGui::End();
}