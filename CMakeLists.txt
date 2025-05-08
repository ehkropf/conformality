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

#add_executable(fftcheck
#    oned_fft.cpp
#)
#
#set_target_properties(fftcheck
#    PROPERTIES INCLUDE_DIRECTORIES "${INCLUDE_DIRS}"
#)

#add_executable(gtest                                                                       target_link_libraries(fftcheck PRIVATE ${FFTW_LIB})
#     main.mm
#     external/imgui/imgui_demo.cpp
#)
#add_executable(pdemo
#    pltdemo.mm
#    external/imgui/imgui_demo.cpp
#    external/implot/implot_demo.cpp
#)
#
## Add Objective-C++ flags where needed.
#set_source_files_properties(main.mm pltdemo.mm external/imgui/backends/imgui_impl_metal.mm
#    PROPERTIES COMPILE_FLAGS "-ObjC++ -fobjc-weak -fobjc-arc"
#)
#set_target_properties(imobjects gtest pdemo
#    PROPERTIES INCLUDE_DIRECTORIES "${IM_INCLUDE_DIRS}"
#)
#target_link_libraries(gtest PRIVATE ${IM_LINK_LIBS} imobjects)
#target_link_libraries(pdemo PRIVATE ${IM_LINK_LIBS} imobjects)
