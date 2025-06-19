#pragma once

#include <memory>

struct GLFWwindow;
class MainWindow;
class Renderer;

class Application
{
private:
    GLFWwindow* mp_window;
    std::unique_ptr<MainWindow> mp_mainWindow;
    std::unique_ptr<Renderer> mp_renderer;
    
    bool m_shouldClose;
    bool m_imguiInitialized;
    
public:
    Application();
    ~Application();
    
    bool initialize();
    void run();
    void shutdown();
    
    void requestClose() { m_shouldClose = true; }
    
private:
    bool setupImGui();
    bool setupWindow();
    void handleEvents();
    void render();
    void cleanup();
    
    static void glfwErrorCallback(int error, const char* description);
};