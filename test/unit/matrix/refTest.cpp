#include <matrixTest.hpp>


TEST_CASE("Matrix ref method", "[ref]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    matrix3.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9})
           .addRow({10, 11, 12});

    SECTION("Normal case with square matrix") {
        Matrix<double> refMatrix = matrix1.ref(tolerance);
        REQUIRE(checkVecApprox(refMatrix.getRow(0), {7, 8, 9}));
        REQUIRE(checkVecApprox(refMatrix.getRow(1), {0, 0.8571428571, 1.7142857143}));
        REQUIRE(checkVecApprox(refMatrix.getRow(2), {0, 0, 0}));
    }

    SECTION("Normal case with rectangular matrix") {
        Matrix<double> refMatrix = matrix2.ref(tolerance);
        REQUIRE(checkVecApprox(refMatrix.getRow(0), {9, 10, 11, 12}));
    }

    SECTION("Normal case with tall matrix") {
        Matrix<double> refMatrix = matrix3.ref(tolerance);
        REQUIRE(checkVecApprox(refMatrix.getRow(0), {10, 11, 12}));
    }
}