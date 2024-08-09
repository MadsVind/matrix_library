# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add Catch2
include(AddCatch2)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Get a list of all test source files
file(GLOB_RECURSE TEST_SOURCES "*.cpp" "../src/*.cpp")

# Add test executable
add_executable(cnn_test ${TEST_SOURCES})

# Set the include directories for the tests
target_include_directories(cnn_test PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Link the required libraries to the test executable
target_link_libraries(cnn_test PRIVATE cnn_lib Catch2::Catch2WithMain)

# Include the CTest module
include(CTest)
include(Catch)

# Discover the tests
catch_discover_tests(cnn_test)