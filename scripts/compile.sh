#!/bin/bash

# Default build type to Release
BUILD_TYPE=Release
PRINT_BENCHMARKS=false
USE_EMSCRIPTEN=false

# Parse command line arguments
while getopts "dbe" opt; do
  case $opt in
    d)
      BUILD_TYPE=Debug
      ;;
    b)
      PRINT_BENCHMARKS=true
      ;;
    e)
      USE_EMSCRIPTEN=true
      ;;
    *)
      echo "Usage: $0 [-d] [-b] [-e]"
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
if [ "$USE_EMSCRIPTEN" == true ]; then
    source ../../emsdk/emsdk_env.sh
    emcmake cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ../..
else
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ../..
fi

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

# Ensure the logs directory exists
mkdir -p ../../logs
echo "Logs directory ensured."

# Run benchmarks if -b flag is used
if [ "$PRINT_BENCHMARKS" == true ]; then
    if [ -f ./bin/cnn_benchmark ]; then
        echo "Running benchmarks directing to ../../logs/benchmark_log.txt ..."
        ./bin/cnn_benchmark --benchmark-samples 100 --out ../../logs/benchmark_log.txt
        if [ -f ../../logs/benchmark_log.txt ]; then
            echo "Benchmark log created successfully at ../../logs/benchmark_log.txt"
        else
            echo "Error: Benchmark log not created."
        fi
    else
        echo "Error: ./bin/cnn_benchmark not found."
    fi
fi