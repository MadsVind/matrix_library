#include <matrixTest.hpp>

TEST_CASE("Matrix Power Operation", "[benchmark][pow]") {
    Matrix<double> matrix1;

    matrix1.addRow({1, 2, 3})
	       .addRow({4, 5, 6})
	       .addRow({7, 8, 9});
    
    BENCHMARK("Calculate Power") {
        return matrix1.pow(sampleSize);
    };
} 