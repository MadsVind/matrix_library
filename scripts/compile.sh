#!/bin/bash

# Default build type to Release
BUILD_TYPE=Release
PRINT_BENCHMARKS=false

# Parse command line arguments
while getopts "db" opt; do
  case $opt in
    d)
      BUILD_TYPE=Debug
      ;;
    b)
      PRINT_BENCHMARKS=true
      ;;
    *)
      echo "Usage: $0 [-d] [-b]"
      exit 1
      ;;
  esac
done

# Set build directory based on build type
if [ "$BUILD_TYPE" == "Debug" ]; then
    BUILD_DIR="build/debug"
else
    BUILD_DIR="build/release"
fi

# Create and navigate to the build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Run CMake to generate the Makefile with the specified build type
# Point CMake to the root directory where CMakeLists.txt is located
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ../..


# Run make to build the project
make

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful."
else
    echo "Compilation failed."
    exit 1
fi

# Conditional execution based on build type
if [ "$BUILD_TYPE" == "Debug" ]; then
    echo "Debug build complete."
    # Debug-specific commands here
else
    # Run tests only in Release mode
    if [ -f ./bin/cnn_test ]; then
        ./bin/cnn_test
        echo "Release build complete."
    else
        echo "Error: ./bin/cnn_test not found."
    fi
fi

# Run the main program
if [ -f ./bin/cnn ]; then
    echo "Running the program..."
    ./bin/cnn
else
    echo "Error: ./bin/cnn not found."
fi

# Run benchmarks if -b flag is used
if [ "$PRINT_BENCHMARKS" == true ]; then
    if [ -f ./bin/cnn_benchmark ]; then
        echo "Running benchmarks..."
        ./bin/cnn_benchmark
    else
        echo "Error: ./bin/cnn_benchmark not found."
    fi
fi