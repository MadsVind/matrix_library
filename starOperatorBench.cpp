#include <matrixTest.hpp>

/*
    TODO: make benchmark and find break point where thread implementation more effective
*/


TEST_CASE("Matrix Multiplication Operator", "[benchmark][operator*]") {
    for (int i = 0; i < sampleSize; ++i) {
        Matrix<double> matrix1 = createIncrementingSquareMatrix<double>(i + 1);
        Matrix<double> matrix2 = matrix1.transpose();

        std::string matrixSizeString = std::to_string(i + 1);
        BENCHMARK("* Operator " + matrixSizeString + " x " + matrixSizeString + "Matrix") {
            return matrix1 * matrix2; 
        };
    }
}