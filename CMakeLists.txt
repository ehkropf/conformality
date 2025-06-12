cmake_minimum_required(VERSION 3.14)
project("Conformality")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wformat")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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
    src/FFTWWrapper.cpp
)

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
    ${CMAKE_SOURCE_DIR}/src/ConformalMap.cpp
    ${CMAKE_SOURCE_DIR}/src/ConformalMapMethod.cpp
    ${CMAKE_SOURCE_DIR}/src/TheodorsenMethod.cpp
)
add_library(util_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/Grid.cpp
    ${CMAKE_SOURCE_DIR}/src/RootFinder.cpp
    ${CMAKE_SOURCE_DIR}/src/StatusManager.cpp
    ${CMAKE_SOURCE_DIR}/src/Types.cpp
)
add_library(domain_objects OBJECT
    ${CMAKE_SOURCE_DIR}/src/Domain.cpp
    ${CMAKE_SOURCE_DIR}/src/BoundaryComponent.cpp
    ${CMAKE_SOURCE_DIR}/src/Boundary.cpp
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
    tests/theodorsen_method_tests.cpp
    $<TARGET_OBJECTS:map_objects>
    tests/grids.cpp
    $<TARGET_OBJECTS:util_objects>
    $<TARGET_OBJECTS:domain_objects>
    tests/root_finder_tests.cpp
    tests/status_manager_tests.cpp
    tests/types_tests.cpp
)
target_link_libraries(test_exec
    GTest::gtest_main
    ${FFTW_LIB}
)

include(GoogleTest)
gtest_discover_tests(test_exec)
