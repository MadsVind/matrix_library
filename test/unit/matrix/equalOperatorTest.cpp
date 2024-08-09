#include <matrixTest.hpp>

TEST_CASE("Does matrix operator== return expected boolean value", "[unit_test]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<double> matrix4;
    Matrix<double> matrix5;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix3.addRow({9, 8, 7})
           .addRow({6, 5, 4})
           .addRow({3, 2, 1});

    // Edge case: Matrix with different dimensions
    matrix4.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    // Edge case: Empty matrix
    matrix5 = Matrix<double>();

    REQUIRE((matrix1 == matrix2) == true);
    REQUIRE((matrix1 == matrix3) == false);
    REQUIRE((matrix1 == matrix4) == false);
    REQUIRE((matrix1 == matrix5) == false);
    REQUIRE((matrix5 == matrix5) == true);
}