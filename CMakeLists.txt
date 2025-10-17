cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

# Project name
project(l4project)



# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the RealSense SDK
find_package(realsense2 REQUIRED)
find_package(glfw3 3.4 REQUIRED)
find_package(OpenGL REQUIRED)
set(OpenCV_DIR /usr/local/opt/opencv/lib/cmake/opencv4)
find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED)
find_package(Freetype REQUIRED)


set(NetCDFCxx_INCLUDE_DIR /usr/local/Cellar/netcdf-cxx/4.3.1_3/include)
set(NetCDFCxx_LIBRARY /usr/local/Cellar/netcdf-cxx/4.3.1_3/lib/libnetcdf-cxx4.1.1.0.dylib)


# Include directories 
include_directories(${realsense2_INCLUDE_DIR})
include_directories(${GLFW_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${NetCDFCxx_INCLUDE_DIR})
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${FREETYPE_INCLUDE_DIRS})


# Link directories
link_directories(${GLFW_LIBRARY_DIRS})
link_directories(${FREETYPE_LIBRARY_DIRS})

# Define source files
set(SOURCES src/main.cpp src/glad.cpp)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE include)

# Link libraries
target_link_libraries(${PROJECT_NAME} ${realsense2_LIBRARY} glfw OpenGL::GL ${OpenCV_LIBS} ${NetCDFCxx_LIBRARY} ${VTK_LIBRARIES} ${FREETYPE_LIBRARIES})

# Optional: Add compiler flags (if needed)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wno-deprecated-declarations)
endif()

# Optional: Install rule
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
