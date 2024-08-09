#include <matrixTest.hpp>

TEST_CASE("Matrix power function works correctly", "[unit_test]") {
    Matrix<double> matrix;
    Matrix<double> expectedMatrix;

    SECTION("Positive integer exponent") {
        // Initialize matrix (2x2)
        matrix.addRow({2, 0})
              .addRow({0, 2});

        // Expected result of matrix^2
        expectedMatrix.addRow({4, 0})
                      .addRow({0, 4});

        REQUIRE(matrix.pow(2) == expectedMatrix);
    }

    SECTION("Negative integer exponent") {
        // Initialize matrix (2x2)
        matrix.addRow({2, 0})
              .addRow({0, 2});

        // Expected result of matrix^-1
        expectedMatrix.addRow({0.5, 0})
                      .addRow({0, 0.5});

        REQUIRE(matrix.pow(-1) == expectedMatrix);
    }

    SECTION("Zero exponent") {
        // Initialize matrix (2x2)
        matrix.addRow({2, 0})
              .addRow({0, 2});

        // Expected result of matrix^0 (identity matrix)
        expectedMatrix.addRow({1, 0})
                      .addRow({0, 1});

        REQUIRE(matrix.pow(0) == expectedMatrix);
    }

    SECTION("Non-square matrix") {
        // Initialize non-square matrix (2x3)
        matrix.addRow({1, 2, 3})
              .addRow({4, 5, 6});

        REQUIRE_THROWS(matrix.pow(2));
    }

    SECTION("Matrix with non-invertible matrix for negative exponent") {
        // Initialize singular matrix (2x2)
        matrix.addRow({1, 2})
              .addRow({2, 4});

        REQUIRE_THROWS(matrix.pow(-1));
    }

    SECTION("Fractional exponent") {
        // Initialize matrix (2x2)
        matrix.addRow({4, 0})
              .addRow({0, 4});

        // Expected result of matrix^0.5 (square root of matrix)
        expectedMatrix.addRow({2, 0})
                      .addRow({0, 2});

        REQUIRE(matrix.pow(0.5) == expectedMatrix);
    }
}