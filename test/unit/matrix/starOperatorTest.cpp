#include <matrixTest.hpp>

TEST_CASE("Matrix multiplication operator* works correctly", "[operator*]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> expectedMatrix;
    std::vector<double> vec1;
    std::vector<double> expectedVec;

    // Initialize matrix1 (3x3)
    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    // Initialize matrix2 (3x3)
    matrix2.addRow({9, 8, 7})
           .addRow({6, 5, 4})
           .addRow({3, 2, 1});

    // Expected result of matrix1 * matrix2
    expectedMatrix.addRow({30, 24, 18})
                  .addRow({84, 69, 54})
                  .addRow({138, 114, 90});

    SECTION("Standard case") {
        // Test matrix multiplication
        REQUIRE((matrix1 * matrix2) == expectedMatrix);
    }
    
    SECTION("Standard matrix vector case") {
        // Initialize vector (3 elements)
        vec1 = {1, 2, 3};

        // Expected result of matrix1 * vec1
        expectedVec = {14, 32, 50};

        // Test matrix-vector multiplication
        REQUIRE((matrix1 * vec1) == expectedVec);
    }

    SECTION("Standard scalar case") {
        // Test scalar multiplication
        double scalar = 2.0;
        Matrix<double> expectedScalarMatrix;
        expectedScalarMatrix.addRow({2, 4, 6})
                            .addRow({8, 10, 12})
                            .addRow({14, 16, 18});

        REQUIRE((matrix1 * scalar) == expectedScalarMatrix);
    }

    SECTION("Matrix A col and matrix B rows doesn't match") { 
        // Edge case: Matrix with different dimensions
        Matrix<double> matrix3;
        matrix3.addRow({1, 2})
               .addRow({3, 4});

        REQUIRE_THROWS(matrix1 * matrix3);
    }

    SECTION("Matrix A col and vec B doesn't match") { 
        // Edge case: Matrix with different dimensions
        std::vector<double> vec2;
        vec2 = {1, 2, 3, 4};

        REQUIRE_THROWS(matrix1 * vec2);
    }
}