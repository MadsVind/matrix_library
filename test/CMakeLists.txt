include(AddCatch2)

file(GLOB_RECURSE TEST_SOURCES "*.cpp" "../src/*.cpp")
add_executable(cnn_tests ${TEST_SOURCES})

# Set the include directories for the tests
target_include_directories(cnn_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Link the required libraries to the test executable
target_link_libraries(cnn_tests PRIVATE cnn_lib Catch2::Catch2WithMain)

# Include the CTest module
include(CTest)
include(Catch)

# Discover the tests
catch_discover_tests(cnn_tests)