#include <matrixTest.hpp>

TEST_CASE("Does matrix transpose return expected matrix", "[unit_test]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<double> matrix4;
    Matrix<double> zeroMatrix(3, 3, 0);

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 4, 7})
           .addRow({2, 5, 8})
           .addRow({3, 6, 9});

    // Edge case: Matrix with different dimensions
    matrix3.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8});

    matrix4.addRow({1, 5})
           .addRow({2, 6})
           .addRow({3, 7})
           .addRow({4, 8});

    zeroMatrix = Matrix<double>(); // Zero matrix

    REQUIRE(matrix1.transpose() == matrix2);
    REQUIRE(matrix3.transpose() == matrix4);
    REQUIRE(matrix4.transpose() == matrix3);
    REQUIRE(zeroMatrix.transpose() == zeroMatrix); // Transpose of zero matrix should still be a zero matrix
}