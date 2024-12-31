#include <matrixTest.hpp>

TEST_CASE("Matrix Row Echelon Form", "[benchmark][ref]") {
    Matrix<double> matrix1;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    BENCHMARK("ref normal matrix") {
        return matrix1.ref();
    };

}