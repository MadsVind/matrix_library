#include <matrixTest.hpp>

TEST_CASE("Does matrix operator+ return expected matrix", "[unit_test]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<double> matrix4;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix3.addRow({2, 4, 6})
           .addRow({8, 10, 12})
           .addRow({14, 16, 18});

    // Edge case: Matrix with different dimensions
    matrix4.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    REQUIRE((matrix1 + matrix2) == matrix3);
    REQUIRE_THROWS_AS(matrix1 + matrix4, std::invalid_argument);
}