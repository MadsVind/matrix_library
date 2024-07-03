#!/bin/bash

# Create a build directory and navigate into it
mkdir -p build
cd build

# Run CMake to generate the Makefile
cmake ..

# Run make to build the project
make

# Run the tests
#ctest

# If the tests were successful, run the program
if [ $? -eq 0 ]; then
    ./cnn
fi