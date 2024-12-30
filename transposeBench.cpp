#include <matrixTest.hpp>

TEST_CASE("Matrix Transpose Operation", "[benchmark][transpose]") {
    Matrix<double> matrix1;

    // Initialize matrix1 (3x3)
    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    BENCHMARK("Tranpose Normal Matrix") {
        return matrix1.transpose();
    };
}