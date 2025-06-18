#include "Renderer.h"

#ifdef __APPLE__
#include "MetalRenderer.h"
#endif

std::unique_ptr<Renderer> Renderer::create()
{
#ifdef __APPLE__
    return std::make_unique<MetalRenderer>();
#else
    // Future: Add OpenGL renderer for other platforms
    return nullptr;
#endif
}