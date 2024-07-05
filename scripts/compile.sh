#!/bin/bash

# Default build type to Release
BUILD_TYPE=Release

# Check if the first argument is -d, if so, set BUILD_TYPE to Debug
if [ "$1" == "-d" ]; then
    BUILD_TYPE=Debug
fi

# Create a build directory and navigate into it
mkdir -p build
cd build

# Run CMake to generate the Makefile with the specified build type
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ..

# Run make to build the project
make

# Conditional execution based on build type
if [ "$BUILD_TYPE" == "Debug" ]; then
    echo "Debug build complete."
    # Debug-specific commands here
else
    # Run tests only in Release mode
    ctest
    echo "Release build complete."
    # Release-specific commands here
fi

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program..."
    # Run the program
    ./cnn
else
    echo "Compilation failed."
fi