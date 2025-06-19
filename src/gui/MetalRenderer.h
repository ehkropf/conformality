#pragma once

#include "Renderer.h"
#include <memory>

class MetalRenderer : public Renderer
{
private:
    void* mp_device;
    void* mp_commandQueue;
    void* mp_layer;
    void* mp_renderPassDescriptor;
    void* mp_currentDrawable;
    GLFWwindow* mp_currentWindow;
    bool m_isShutdown;
    
public:
    MetalRenderer();
    ~MetalRenderer() override;
    
    bool initialize(GLFWwindow* window) override;
    void beginFrame() override;
    void endFrame() override;
    void shutdown() override;
    
private:
    bool setupMetal(GLFWwindow* window);
};