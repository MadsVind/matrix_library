#include <matrixTest.hpp>

#define MATRIX_UTIL

const size_t sampleSize = 1000;

TEST_CASE("Does isSquare() return the expected boolean value", "[isSquare]") {
    for (int i = 1; i < sampleSize; ++i) {
        int j = sampleSize - i;
        Matrix<double> probNonSquare = Matrix<double>(i, j, 0);
        REQUIRE(probNonSquare.isSquare() == (i == j)); 
    }

    for (int i = 1; i < sampleSize; ++i) {
        Matrix<double> square = Matrix<double>(i, i, 0);
        REQUIRE(square.isSquare() == true); 
    }
}

TEST_CASE("does isSameSize() return the expected boolean value", "[isSameSize]") {
    for (int i = 1; i < sampleSize; ++i) {
        int j = sampleSize - i;
        Matrix<double> ij = Matrix<double>(i, j, 0);
        Matrix<double> ji = Matrix<double>(j, i, 0);
        REQUIRE(ij.isSameSize(ji) == (i == j)); 
    }

    for (int i = 1; i < sampleSize; ++i) {
        int j = sampleSize - i;
        Matrix<double> ij = Matrix<double>(i, j, 0);
        Matrix<double> ij2 = Matrix<double>(i, j, 0);
        REQUIRE(ij.isSameSize(ij2) == true); 
    }
}