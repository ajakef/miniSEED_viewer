cmake_minimum_required(VERSION 3.16)
project(miniseed_viewer LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(TP ${CMAKE_SOURCE_DIR}/third_party)

find_package(Threads REQUIRED)

# --- libmseed (reads miniSEED) ---
file(GLOB MSEED_SRC ${TP}/libmseed/*.c)
add_library(mseed STATIC ${MSEED_SRC})
target_include_directories(mseed PUBLIC ${TP}/libmseed)
target_compile_options(mseed PRIVATE -w)

# --- the part of the program that has no graphics ---
add_library(core STATIC src/data.cpp src/dsp.cpp src/stats.cpp src/settings.cpp)
target_include_directories(core PUBLIC ${TP}/libmseed ${TP}/pocketfft)
target_link_libraries(core PUBLIC mseed Threads::Threads)
if(NOT MSVC)
  target_compile_options(core PRIVATE -Wall -Wextra)
endif()

add_executable(test_core tests/test_core.cpp)
target_link_libraries(test_core core)

# --- the program ---
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "" FORCE)  # X11 only: also works on Wayland desktops through XWayland
add_subdirectory(${TP}/glfw)
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)

add_library(imgui STATIC
  ${TP}/imgui/imgui.cpp ${TP}/imgui/imgui_draw.cpp ${TP}/imgui/imgui_tables.cpp ${TP}/imgui/imgui_widgets.cpp
  ${TP}/imgui/imgui_impl_glfw.cpp ${TP}/imgui/imgui_impl_opengl3.cpp
  ${TP}/implot/implot.cpp ${TP}/implot/implot_items.cpp)
target_include_directories(imgui PUBLIC ${TP}/imgui ${TP}/implot)
target_link_libraries(imgui PUBLIC glfw OpenGL::GL)
target_compile_options(imgui PRIVATE -w)

add_executable(miniseed_viewer src/main.cpp src/app.cpp src/app_views.cpp src/app_specgram.cpp)
target_link_libraries(miniseed_viewer core imgui Threads::Threads)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_options(miniseed_viewer PRIVATE -static-libstdc++ -static-libgcc)
endif()
if(NOT MSVC)
  target_compile_options(miniseed_viewer PRIVATE -Wall -Wextra)
endif()
