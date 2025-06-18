#include "Application.h"
#include "MainWindow.h"
#include "Renderer.h"
#include "imgui.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <stdio.h>

Application::Application()
    : mp_window{nullptr}
    , mp_mainWindow{nullptr}
    , mp_renderer{nullptr}
    , m_shouldClose{false}
{
}

Application::~Application()
{
    shutdown();
}

bool Application::initialize()
{
    if (!setupWindow())
    {
        return false;
    }
    
    if (!setupImGui())
    {
        return false;
    }
    
    mp_renderer = Renderer::create();
    if (!mp_renderer || !mp_renderer->initialize(mp_window))
    {
        return false;
    }
    
    mp_mainWindow = std::make_unique<MainWindow>();
    if (!mp_mainWindow->initialize())
    {
        return false;
    }
    
    return true;
}

void Application::run()
{
    while (!glfwWindowShouldClose(mp_window) && !m_shouldClose)
    {
        handleEvents();
        render();
    }
}

void Application::shutdown()
{
    if (mp_mainWindow)
    {
        mp_mainWindow->shutdown();
        mp_mainWindow.reset();
    }
    
    if (mp_renderer)
    {
        mp_renderer->shutdown();
        mp_renderer.reset();
    }
    
    cleanup();
}

bool Application::setupImGui()
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
    
    ImGui::StyleColorsDark();
    
    // Use default font for now
    io.Fonts->AddFontDefault();
    
    return true;
}

bool Application::setupWindow()
{
    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit())
    {
        return false;
    }
    
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    mp_window = glfwCreateWindow(1280, 720, "Conformality Mapping Tool", nullptr, nullptr);
    if (!mp_window)
    {
        return false;
    }
    
    return true;
}

void Application::handleEvents()
{
    glfwPollEvents();
    
    // Update layer drawable size for Metal
    if (mp_renderer)
    {
        int width, height;
        glfwGetFramebufferSize(mp_window, &width, &height);
        // The renderer handles platform-specific details
    }
}

void Application::render()
{
    if (mp_renderer)
    {
        mp_renderer->beginFrame();
    }
    
    if (mp_mainWindow)
    {
        mp_mainWindow->render();
    }
    
    if (mp_renderer)
    {
        mp_renderer->endFrame();
    }
}

void Application::cleanup()
{
    ImGui::DestroyContext();
    
    if (mp_window)
    {
        glfwDestroyWindow(mp_window);
        mp_window = nullptr;
    }
    
    glfwTerminate();
}

void Application::glfwErrorCallback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}