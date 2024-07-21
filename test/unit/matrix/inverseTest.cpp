#include <matrixTest.hpp>

TEST_CASE("Matrix inverse method", "[inverse]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<double> matrix4;
    Matrix<double> matrix5;

    matrix1.addRow({2, 0, 2, 0.6})
           .addRow({3, 3, 4, -2})
           .addRow({5, 5, 4, 2})
           .addRow({-1, -2, 3.4, -1});

    matrix2.addRow({1, 0, 0})
           .addRow({0, 1, 0})
           .addRow({0, 0, 1});

    matrix3.addRow({0, 2, 2})
           .addRow({1, 2, 3})
           .addRow({4, 4, 4});

    matrix4.addRow({0, 2, 2})
           .addRow({0, 2, 3})
           .addRow({0, 4, 4});

    matrix5.addRow({2, 2, 2, 3})
           .addRow({3, 2, 3, 2})
           .addRow({1, 4, 4, 1});

    SECTION("Simple matrix") {
        Matrix<double> expectedInverse;
        expectedInverse.addRow({0.653333, 0.141667, -0.207667, -0.306667})
                    .addRow({-0.52, 0.025, 0.201, 0.04})
                    .addRow({-0.133333, -0.0416667, 0.131667, 0.266667})
                    .addRow({-0.0666667, -0.333333, 0.253333, 0.133333});

        REQUIRE(checkMatrixApprox(matrix1.inverse(), expectedInverse));
    }
    SECTION("Identity matrix") {
        Matrix<double> expectedInverse;
        expectedInverse.addRow({1, 0, 0})
                       .addRow({0, 1, 0})
                       .addRow({0, 0, 1});

        REQUIRE(checkMatrixApprox(matrix2.inverse(), expectedInverse));
    }

    SECTION("Matrix requiring row swaps") {
        Matrix<double> expectedInverse;
        expectedInverse.addRow({-0.5, 0, 0.25})
                       .addRow({1, -1, 0.25})
                       .addRow({-0.5, 1, -0.25});

        REQUIRE(checkMatrixApprox(matrix3.inverse(), expectedInverse));
    }

    SECTION("Matrix with empty column") {
        REQUIRE_THROWS(matrix4.inverse());
    }

    SECTION("Nonsquare matrix") {
        REQUIRE_THROWS(matrix4.inverse());
    }
}