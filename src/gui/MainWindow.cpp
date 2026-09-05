#include "MainWindow.h"
#include "VisualizationPanel.h"
#include "GuiController.h"
#include "Application.h"
#include "../core/MethodInfo.h"
#include "imgui.h"
#include "implot.h"

using conformality::formatMethodInfoValue;

MainWindow::MainWindow()
    : m_showDemoWindow{false}
    , m_showAboutDialog{false}
    , mp_visualizationPanel{nullptr}
    , mp_controller{nullptr}
    , m_statusLog{"Ready"}
    , mp_application{nullptr}
    , m_showGrid{true}
    , m_gridDensity{8}
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

    return true;
}

void MainWindow::render()
{
    if (mp_controller)
    {
        mp_controller->update();
    }

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
            bool computing = mp_controller && mp_controller->isComputing();
            if (computing) ImGui::BeginDisabled();
            if (ImGui::MenuItem("New"))
            {
                if (mp_controller)
                {
                    mp_controller->reset();
                }
            }
            if (computing) ImGui::EndDisabled();
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

        if (ImGui::BeginMenu("Examples"))
        {
            bool computing = mp_controller && mp_controller->isComputing();
            if (computing) ImGui::BeginDisabled();

            if (ImGui::MenuItem("Thesis Ex 3 (Identity, m=4)"))
            {
                if (mp_controller) mp_controller->loadThesisExample(3);
            }
            if (ImGui::MenuItem("Thesis Ex 5 (Ellipses, m=3)"))
            {
                if (mp_controller) mp_controller->loadThesisExample(5);
            }
            if (ImGui::MenuItem("Thesis Ex 2 (Mixed, m=4)"))
            {
                if (mp_controller) mp_controller->loadThesisExample(2);
            }
            if (ImGui::MenuItem("Thesis Ex 4 (High connectivity, m=7)"))
            {
                if (mp_controller) mp_controller->loadThesisExample(4);
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Reset"))
            {
                if (mp_controller) mp_controller->reset();
            }

            if (computing) ImGui::EndDisabled();
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
        // Snapshot live progress once per frame for consistent display
        if (mp_controller && mp_controller->isComputing())
        {
            mp_controller->getLiveProgress(m_cachedLiveIter, m_cachedLiveResidual);
        }

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

    if (ImGui::BeginChild("ControlPanel", ImVec2(control_width, -m_statusPanelHeight), true))
    {
        ImGui::Text("Control Panel");
        ImGui::Separator();

        // Map info
        if (mp_controller && mp_controller->hasMap())
        {
            ImGui::TextWrapped("%s", mp_controller->getMapDescription().c_str());
        }
        else
        {
            ImGui::TextDisabled("No map loaded");
            ImGui::TextDisabled("Use Examples menu to load a preset");
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

        if (ImGui::SliderInt("Grid Density", &m_gridDensity, 4, 16))
        {
            if (mp_visualizationPanel)
            {
                mp_visualizationPanel->setGridDensity(m_gridDensity);
            }
        }

        ImGui::Separator();

        // Computation results from MethodInfo
        if (mp_controller)
        {
            const auto& info = mp_controller->getLastMethodInfo();
            if (!info.results.empty())
            {
                for (const auto& field : info.results)
                {
                    auto display = formatMethodInfoValue(field.value);
                    ImGui::Text("%s: %s", field.label.c_str(), display.c_str());
                }
                ImGui::Separator();
            }
        }

        // Compute button
        bool isComputing = mp_controller ? mp_controller->isComputing() : false;
        bool hasMap = mp_controller ? mp_controller->hasMap() : false;

        if (isComputing)
        {
            if (m_cachedLiveIter > 0)
            {
                ImGui::Text("Iteration: %d", m_cachedLiveIter);
                ImGui::Text("Residual: %.2e", m_cachedLiveResidual);
            }
            else
            {
                ImGui::Text("%s", m_computationPhase.c_str());
            }

            if (ImGui::Button("Cancel", ImVec2(-1, 0)))
            {
                if (mp_controller)
                {
                    mp_controller->cancelComputation();
                }
            }
        }
        else
        {
            if (!m_computationPhase.empty())
            {
                ImGui::Text("%s", m_computationPhase.c_str());
            }

            if (!hasMap) ImGui::BeginDisabled();
            if (ImGui::Button("Compute Mapping", ImVec2(-1, 0)))
            {
                if (mp_controller)
                {
                    if (mp_controller->computeMapping())
                    {
                        m_computationPhase = "Computing...";
                    }
                }
            }
            if (!hasMap) ImGui::EndDisabled();
        }
    }
    ImGui::EndChild();
}

void MainWindow::renderVisualizationPanel()
{
    if (ImGui::BeginChild("VisualizationPanel", ImVec2(0, -m_statusPanelHeight), true))
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
    // ImGuiChildFlags_ResizeY lets the user drag the panel's top border to resize it. BeginChild
    // manages the size internally once resized, so the requested height below only matters on the
    // first frame; the actual current height is read back via GetWindowSize() after EndChild so
    // sibling panels (ControlPanel, VisualizationPanel) can reserve the correct remaining space.
    if (ImGui::BeginChild("StatusPanel", ImVec2(0, m_statusPanelHeight),
                           ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeY))
    {
        // Info bar: live progress, errors, results, performance
        if (mp_controller)
        {
            if (mp_controller->isComputing())
            {
                if (m_cachedLiveIter > 0)
                {
                    ImGui::Text("Iter: %d | Residual: %.2e", m_cachedLiveIter, m_cachedLiveResidual);
                }
                else
                {
                    ImGui::Text("Computing...");
                }
            }
            else
            {
                if (!mp_controller->getLastErrorMessage().empty())
                {
                    ImGui::Text("Error: %s", mp_controller->getLastErrorMessage().c_str());
                }
                else
                {
                    const auto& info = mp_controller->getLastMethodInfo();
                    for (const auto& field : info.results)
                    {
                        auto display = formatMethodInfoValue(field.value);
                        ImGui::Text("%s: %s  |", field.label.c_str(), display.c_str());
                        ImGui::SameLine();
                    }
                    ImGui::Text("Performance: %.3f ms/frame (%.1f FPS)",
                               1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
                }
            }
        }
        else
        {
            ImGui::Text("Performance: %.3f ms/frame (%.1f FPS)",
                       1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        }

        ImGui::Separator();

        // Scrolling status log
        if (ImGui::BeginChild("StatusLog", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar))
        {
            for (const auto& entry : m_statusLog)
            {
                ImGui::TextUnformatted(entry.c_str());
            }

            if (m_statusLogScrollToBottom)
            {
                ImGui::SetScrollHereY(1.0f);
                m_statusLogScrollToBottom = false;
            }
        }
        ImGui::EndChild();

        // Persist the current (possibly user-resized) height for next frame's sibling layout.
        m_statusPanelHeight = ImGui::GetWindowSize().y;
    }
    ImGui::EndChild();
}

void MainWindow::onStatusUpdate(const std::string& message)
{
    // Cap log size to avoid unbounded growth
    static constexpr size_t kMaxLogEntries = 1000;
    if (m_statusLog.size() >= kMaxLogEntries)
    {
        m_statusLog.erase(m_statusLog.begin(), m_statusLog.begin() + static_cast<long>(m_statusLog.size() / 2));
    }

    m_statusLog.push_back(message);
    m_statusLogScrollToBottom = true;
    m_computationPhase = message;
}

void MainWindow::onComputationComplete()
{
    if (!mp_controller)
    {
        return;
    }

    if (mp_controller->wasCancelled())
    {
        m_computationPhase = "Cancelled";
    }
    else if (mp_controller->wasLastComputationSuccessful())
    {
        const auto& info = mp_controller->getLastMethodInfo();
        bool converged = false;
        int iterations = 0;
        for (const auto& field : info.results)
        {
            if (field.label == "Converged")
            {
                if (auto* p = std::get_if<bool>(&field.value))
                    converged = *p;
            }
            else if (field.label == "Iterations")
            {
                if (auto* p = std::get_if<int>(&field.value))
                    iterations = *p;
            }
        }

        if (converged)
        {
            m_computationPhase = "Converged (" + std::to_string(iterations) + " iterations)";
        }
        else
        {
            m_computationPhase = "Completed (" + std::to_string(iterations) + " iterations, not converged)";
        }
    }
    else
    {
        m_computationPhase = "Failed: " + mp_controller->getLastErrorMessage();
    }
}

void MainWindow::renderAboutDialog()
{
    if (m_showAboutDialog)
    {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f),
                                ImGuiCond_Always, ImVec2(0.5f, 0.5f));

        if (ImGui::Begin("About Conformality Mapping Tool", &m_showAboutDialog,
                         ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize))
        {
            ImGui::Text("Conformality Mapping Tool");
            ImGui::Text("Version: 1/∞");
            ImGui::Separator();
            ImGui::Text("A tool for visualizing conformal mappings.");
            ImGui::Separator();

            if (ImGui::Button("Close", ImVec2(120, 0)))
            {
                m_showAboutDialog = false;
            }
        }
        ImGui::End();
    }
}
