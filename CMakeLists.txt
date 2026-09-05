cmake_minimum_required(VERSION 3.14)
project("Conformality")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wformat")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#####
# spdlog logging library (added early so object libraries can use it)
set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
set(SPDLOG_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/spdlog)

#####
# Libraries, etc.
set(IM_INCLUDE_DIRS
    /opt/local/include
    ${CMAKE_SOURCE_DIR}/external/imgui
    ${CMAKE_SOURCE_DIR}/external/imgui/backends
    ${CMAKE_SOURCE_DIR}/external/implot
)

find_library(FFTW_LIB fftw3 HINTS /opt/local/lib /usr/local/lib REQUIRED)
if(FFTW_LIB)
    get_filename_component(LIB_DIR "${FFTW_LIB}" DIRECTORY)
    get_filename_component(INSTALL_PREFIX "${LIB_DIR}" DIRECTORY)
    find_path(FFTW_INC_DIR
        NAMES fftw3.h
        PATHS "${INSTALL_PREFIX}/include"
        NO_DEFAULT_PATH
    )
    if(FFTW_INC_DIR)
        include_directories(${FFTW_INC_DIR})
    else()
        message(FATAL_ERROR "Unable to find path for fftw3.h")
    endif()
endif()
add_library(fftw_wrapper OBJECT
    src/numerics/FFTWWrapper.cpp
)

#####
# Eigen linear algebra library
find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Dense
    PATHS /opt/local/include/eigen3 /usr/local/include/eigen3
    NO_DEFAULT_PATH
)
if(EIGEN3_INCLUDE_DIR)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Unable to find Eigen3 installation")
endif()

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

add_library(imobjects OBJECT
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_draw.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_tables.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_widgets.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/imgui_demo.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/backends/imgui_impl_glfw.cpp
    ${CMAKE_SOURCE_DIR}/external/imgui/backends/imgui_impl_metal.mm
    ${CMAKE_SOURCE_DIR}/external/implot/implot.cpp
    ${CMAKE_SOURCE_DIR}/external/implot/implot_items.cpp
)
set_property(TARGET imobjects
    PROPERTY INCLUDE_DIRECTORIES ${IM_INCLUDE_DIRS}
)

#####
# Sources
add_library(map_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/core/ConformalMap.cpp
    ${CMAKE_SOURCE_DIR}/src/methods/ConformalMapMethod.cpp
    ${CMAKE_SOURCE_DIR}/src/methods/FornbergMC.cpp
    ${CMAKE_SOURCE_DIR}/src/methods/FornbergMCConfiguration.cpp
    ${CMAKE_SOURCE_DIR}/src/methods/PMatrixBuilder.cpp
)
target_link_libraries(map_objects PUBLIC spdlog::spdlog_header_only)

add_library(util_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/numerics/Grid.cpp
    ${CMAKE_SOURCE_DIR}/src/numerics/RootFinder.cpp
    ${CMAKE_SOURCE_DIR}/src/numerics/CGSolver.cpp
    ${CMAKE_SOURCE_DIR}/src/core/StatusManager.cpp
    ${CMAKE_SOURCE_DIR}/src/core/Types.cpp
)
target_link_libraries(util_objects PUBLIC spdlog::spdlog_header_only)

add_library(domain_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/domains/Domain.cpp
    ${CMAKE_SOURCE_DIR}/src/domains/BoundaryComponent.cpp
    ${CMAKE_SOURCE_DIR}/src/domains/Boundary.cpp
    ${CMAKE_SOURCE_DIR}/src/domains/FornbergCanonicalDomain.cpp
    ${CMAKE_SOURCE_DIR}/src/domains/InvertedEllipseComponent.cpp
    ${CMAKE_SOURCE_DIR}/src/domains/SplineBoundaryComponent.cpp
)
target_link_libraries(domain_objects PUBLIC spdlog::spdlog_header_only)

add_library(example_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/examples/BoundaryHelpers.cpp
    ${CMAKE_SOURCE_DIR}/src/examples/ThesisExamples.cpp
)
target_link_libraries(example_objects PUBLIC spdlog::spdlog_header_only)

#####
# GUI objects
add_library(gui_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/gui/Application.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/MainWindow.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/Renderer.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/MetalRenderer.mm
    ${CMAKE_SOURCE_DIR}/src/gui/VisualizationPanel.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/GuiController.cpp
)
target_include_directories(gui_objects PRIVATE ${IM_INCLUDE_DIRS})
target_link_libraries(gui_objects PUBLIC spdlog::spdlog_header_only)

# GUI executable
add_executable(conformality_gui
    ${CMAKE_SOURCE_DIR}/src/gui/main.cpp
    $<TARGET_OBJECTS:gui_objects>
    $<TARGET_OBJECTS:imobjects>
    $<TARGET_OBJECTS:map_objects>
    $<TARGET_OBJECTS:util_objects>
    $<TARGET_OBJECTS:domain_objects>
    $<TARGET_OBJECTS:example_objects>
    $<TARGET_OBJECTS:fftw_wrapper>
)
target_link_libraries(conformality_gui
    ${IM_LINK_LIBS}
    ${FFTW_LIB}
    spdlog::spdlog_header_only
)

# CLI executable
add_executable(conformality_cli
    ${CMAKE_SOURCE_DIR}/src/cli/main.cpp
    $<TARGET_OBJECTS:map_objects>
    $<TARGET_OBJECTS:util_objects>
    $<TARGET_OBJECTS:domain_objects>
    $<TARGET_OBJECTS:example_objects>
    $<TARGET_OBJECTS:fftw_wrapper>
)
target_link_libraries(conformality_cli
    ${FFTW_LIB}
    spdlog::spdlog_header_only
)

#####
# Testing section
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/90a41521142c978131f38c6da07b4eb96a9f1ff6.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(test_exec
    tests/complex_rep.cpp
    tests/boundary_components.cpp
    tests/boundary.cpp
    tests/simply_connected.cpp
    tests/starlike_domain.cpp
    tests/elliptical_domain.cpp
    tests/fftw_wrapper.cpp
    $<TARGET_OBJECTS:fftw_wrapper>
    tests/conformal_map_tests.cpp
    tests/fornberg_mc.cpp
    tests/pmatrixbuilder.cpp
    tests/cgsolver.cpp
    tests/polyval.cpp
    $<TARGET_OBJECTS:map_objects>
    tests/grids.cpp
    $<TARGET_OBJECTS:util_objects>
    $<TARGET_OBJECTS:domain_objects>
    $<TARGET_OBJECTS:example_objects>
    tests/root_finder_tests.cpp
    tests/status_manager_tests.cpp
    tests/types_tests.cpp
    tests/constructive_validation.cpp
    tests/inverted_ellipse_tests.cpp
    tests/boundary_helpers_tests.cpp
    tests/thesis_examples_tests.cpp
    tests/spline_boundary_tests.cpp
    tests/method_info_tests.cpp
)
target_link_libraries(test_exec
    GTest::gtest_main
    ${FFTW_LIB}
    spdlog::spdlog_header_only
)

target_compile_definitions(test_exec PRIVATE TESTING)

include(GoogleTest)
gtest_discover_tests(test_exec)
