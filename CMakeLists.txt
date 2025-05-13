cmake_minimum_required(VERSION 3.14)
project("Conformality")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wformat")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(IM_INCLUDE_DIRS
    /opt/local/include
    ${CMAKE_SOURCE_DIR}/external/imgui
    ${CMAKE_SOURCE_DIR}/external/imgui/backends
    ${CMAKE_SOURCE_DIR}/external/implot
)

find_library(FFTW_LIB fftw3 HINTS /opt/local/lib /usr/local/lib REQUIRED)
find_library(GLFW_LIB glfw.3 HINTS /opt/local/lib /usr/local/lib REQUIRED)
find_library(METAL_LIB Metal REQUIRED)
find_library(METALKIT_LIB MetalKit REQUIRED)
find_library(COCOA_LIB Cocoa REQUIRED)
find_library(IOKIT_LIB IOKit REQUIRED)
find_library(COREVID_LIB CoreVideo REQUIRED)
find_library(QUARTZ_LIB QuartzCore REQUIRED)

set(IM_LINK_LIBS
    ${GLFW_LIB}
    ${METAL_LIB}
    ${METALKIT_LIB}
    ${COCOA_LIB}
    ${IOKIT_LIB}
    ${COREVID_LIB}
    ${QUARTZ_LIB}
)

add_library(imobjects
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_draw.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_tables.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_widgets.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/backends/imgui_impl_glfw.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/backends/imgui_impl_metal.mm
    ${CMAKE_SOURCE_DIR}/external/implot/implot.cpp
    ${CMAKE_SOURCE_DIR}/external/implot/implot_items.cpp
)
set_property(TARGET imobjects
    PROPERTY INCLUDE_DIRECTORIES ${IM_INCLUDE_DIRS}
)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/90a41521142c978131f38c6da07b4eb96a9f1ff6.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(test_exec
    tests/complex_rep.cpp
    tests/boundary_components.cpp
    tests/starlike_domain.cpp
    tests/elliptical_domain.cpp
)
target_link_libraries(test_exec GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_exec)
