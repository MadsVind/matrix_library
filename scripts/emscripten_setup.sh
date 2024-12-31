#!/bin/bash

# Update and install dependencies
sudo apt-get update
sudo apt-get install -y git cmake python3 python3-pip

# Install Emscripten SDK
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk

# Fetch the latest version of the Emscripten SDK
./emsdk install latest

# Make the "latest" SDK active for the current user
./emsdk activate latest

# Activate PATH and other environment variables in the current terminal
source ./emsdk_env.sh

# Verify installation
emcc -v

echo "Emscripten setup is complete."