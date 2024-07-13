#include <matrixTest.hpp>

TEST_CASE("Does unary matrix operator- return expected matrix", "[unary operator-]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> zeroMatrix;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({-1, -2, -3})
           .addRow({-4, -5, -6})
           .addRow({-7, -8, -9});

    zeroMatrix = Matrix<double>(3, 3, 0); // Zero matrix

    REQUIRE((-matrix1) == matrix2);
    REQUIRE((-zeroMatrix) == zeroMatrix); // Unary minus on zero matrix should still be a zero matrix
}