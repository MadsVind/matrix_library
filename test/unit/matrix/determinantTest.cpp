#include <matrixTest.hpp>

TEST_CASE("Matrix determinant method", "[determinant]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;

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

    SECTION("Simple matrix") {
        REQUIRE(matrix1.determinant() == -120);
    }

    SECTION("Identity matrix") {
        REQUIRE(matrix2.determinant() == 1);
    }

    SECTION("Matrix requiring row swaps") {  
        REQUIRE(matrix3.determinant() == 8);
    }
}