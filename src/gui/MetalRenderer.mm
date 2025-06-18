#include "MetalRenderer.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_metal.h"

#define GLFW_INCLUDE_NONE
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#import <Metal/Metal.h>
#import <QuartzCore/QuartzCore.h>

MetalRenderer::MetalRenderer()
    : mp_device{nullptr}
    , mp_commandQueue{nullptr}
    , mp_layer{nullptr}
    , mp_renderPassDescriptor{nullptr}
    , mp_currentDrawable{nullptr}
    , mp_currentWindow{nullptr}
{
}

MetalRenderer::~MetalRenderer()
{
    shutdown();
}

bool MetalRenderer::initialize(GLFWwindow* window)
{
    mp_currentWindow = window;
    
    if (!setupMetal(window))
    {
        return false;
    }
    
    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOther(window, true);
    
    id<MTLDevice> device = (__bridge id<MTLDevice>)mp_device;
    ImGui_ImplMetal_Init(device);
    
    return true;
}

void MetalRenderer::beginFrame()
{
    CAMetalLayer* layer = (__bridge CAMetalLayer*)mp_layer;
    MTLRenderPassDescriptor* renderPassDescriptor = (__bridge MTLRenderPassDescriptor*)mp_renderPassDescriptor;
    
    // Get current framebuffer size and update drawable size
    int width, height;
    glfwGetFramebufferSize(mp_currentWindow, &width, &height);
    layer.drawableSize = CGSizeMake(width, height);
    
    // Get the drawable and set up the texture for ImGui
    id<CAMetalDrawable> drawable = [layer nextDrawable];
    renderPassDescriptor.colorAttachments[0].texture = drawable.texture;
    
    // Store drawable for endFrame
    mp_currentDrawable = (void*)CFBridgingRetain(drawable);
    
    ImGui_ImplMetal_NewFrame(renderPassDescriptor);
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void MetalRenderer::endFrame()
{
    ImGui::Render();
    
    id<MTLCommandQueue> commandQueue = (__bridge id<MTLCommandQueue>)mp_commandQueue;
    MTLRenderPassDescriptor* renderPassDescriptor = (__bridge MTLRenderPassDescriptor*)mp_renderPassDescriptor;
    id<CAMetalDrawable> drawable = (__bridge id<CAMetalDrawable>)mp_currentDrawable;
    
    @autoreleasepool
    {
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        
        // Configure render pass descriptor for this frame
        renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(0.45, 0.55, 0.60, 1.0);
        renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
        renderPassDescriptor.colorAttachments[0].storeAction = MTLStoreActionStore;
        
        id<MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
        [renderEncoder pushDebugGroup:@"ImGui Conformality"];
        
        ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), commandBuffer, renderEncoder);
        
        [renderEncoder popDebugGroup];
        [renderEncoder endEncoding];
        
        [commandBuffer presentDrawable:drawable];
        [commandBuffer commit];
    }
    
    // Release the drawable
    if (mp_currentDrawable)
    {
        CFRelease(mp_currentDrawable);
        mp_currentDrawable = nullptr;
    }
}

void MetalRenderer::shutdown()
{
    if (mp_renderPassDescriptor)
    {
        CFRelease(mp_renderPassDescriptor);
        mp_renderPassDescriptor = nullptr;
    }
    
    ImGui_ImplMetal_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    
    if (mp_commandQueue)
    {
        CFRelease(mp_commandQueue);
        mp_commandQueue = nullptr;
    }
    
    if (mp_device)
    {
        CFRelease(mp_device);
        mp_device = nullptr;
    }
    
    if (mp_layer)
    {
        CFRelease(mp_layer);
        mp_layer = nullptr;
    }
}

bool MetalRenderer::setupMetal(GLFWwindow* window)
{
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device)
    {
        return false;
    }
    
    id<MTLCommandQueue> commandQueue = [device newCommandQueue];
    if (!commandQueue)
    {
        return false;
    }
    
    NSWindow* nswin = glfwGetCocoaWindow(window);
    CAMetalLayer* layer = [CAMetalLayer layer];
    layer.device = device;
    layer.pixelFormat = MTLPixelFormatBGRA8Unorm;
    
    // Set initial drawable size
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    layer.drawableSize = CGSizeMake(width, height);
    
    nswin.contentView.layer = layer;
    nswin.contentView.wantsLayer = YES;
    
    MTLRenderPassDescriptor* renderPassDescriptor = [MTLRenderPassDescriptor new];
    
    // Store as void* for cross-platform compatibility
    mp_device = (void*)CFBridgingRetain(device);
    mp_commandQueue = (void*)CFBridgingRetain(commandQueue);
    mp_layer = (void*)CFBridgingRetain(layer);
    mp_renderPassDescriptor = (void*)CFBridgingRetain(renderPassDescriptor);
    
    return true;
}