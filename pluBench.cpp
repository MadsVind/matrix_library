#include <matrixTest.hpp>

// Test cases for decompPLU method
TEST_CASE("PLU Matrix Decompotition", "[benchmark][plu]") {
    Matrix<double> matrix1;

    matrix1.addRow({2, 0, 2, 0.6})
           .addRow({3, 3, 4, -2})
           .addRow({5, 5, 4, 2})
           .addRow({-1, -2, 3.4, -1});

    BENCHMARK("PLU Matrix Decompotition") {
        return matrix1.plu();
    };

}