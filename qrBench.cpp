#include <matrixTest.hpp>

TEST_CASE("Matrix QR Decomposition", "[benchmark][qr]") {
    Matrix<double> matrix1;

    matrix1.addRow({1, 1, 0})
           .addRow({1, 0, 1})
           .addRow({0, 1, 1});

    Matrix<double> sparseMatrix;
    sparseMatrix.addRow({1, 0, 0})
                .addRow({0, 0, 2})
                .addRow({0, 3, 0});

    BENCHMARK("QR Normal Matrix") {
        return matrix1.qr();
    };

    BENCHMARK("QR Sparse Matrix") {
        return sparseMatrix.qr();
    };

}