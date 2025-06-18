#pragma once

#include <memory>

struct GLFWwindow;

class Renderer
{
public:
    virtual ~Renderer() = default;
    
    virtual bool initialize(GLFWwindow* window) = 0;
    virtual void beginFrame() = 0;
    virtual void endFrame() = 0;
    virtual void shutdown() = 0;
    
    static std::unique_ptr<Renderer> create();
};