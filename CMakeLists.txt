# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(cnn VERSION 1.0.0)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Enable testing
enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Include the directories
include_directories("include")
file(GLOB_RECURSE ALL_INCLUDE_DIRS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.hpp" "src/**/*.hpp" "include/**/*.hpp")
foreach(DIR ${ALL_INCLUDE_DIRS})
    get_filename_component(DIR_PATH ${DIR} PATH)
    include_directories(${DIR_PATH})
endforeach()

# Get a list of all source files in the src directory and its subdirectories
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/**/*.cpp" "include/**/*.cpp")

# Add a library that is built from the source files
add_library(cnn_lib ${SOURCES})

# Link the required libraries to the library
target_link_libraries(cnn_lib PUBLIC)

# Check if Emscripten is being used
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".js")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WASM=1 -s MODULARIZE=1 -s EXPORT_NAME=createModule")
endif()

# Add an executable that is built from main.cpp
add_executable(cnn main.cpp)

# Link the library to the executable
target_link_libraries(cnn PRIVATE cnn_lib)

add_subdirectory(test)