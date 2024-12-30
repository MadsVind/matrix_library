#include <matrixTest.hpp>

TEST_CASE("Matrix Minus Operator", "[benchmark][operator-]") {
       Matrix<double> matrix1;
       Matrix<double> matrix2;

       matrix1.addRow({1, 2, 3})
              .addRow({4, 5, 6})
              .addRow({7, 8, 9});

       matrix2.addRow({1, 2, 3})
              .addRow({4, 5, 6})
              .addRow({7, 8, 9});

       BENCHMARK("Matrix Minus Operator") {
              return matrix1 - matrix2;
       };
}