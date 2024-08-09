#include <matrixTest.hpp>

TEST_CASE("Does getIdentity() create identity matrix of apropiate size", "[unit_test]") {
    const size_t matrixTestSize = 9;

    for (int i = 1; i < sampleSize; ++i) {
        Matrix<double> I = Matrix<double>::getIdentity(i);
        REQUIRE(I.getColAmount() == i);
        REQUIRE(I.getRowAmount() == i);
    }
    Matrix<double> I = Matrix<double>::getIdentity(matrixTestSize);
    for (int i = 0; i < matrixTestSize; ++i) {
        for (int j = 0; j < matrixTestSize; ++j) {
            if (i == j) REQUIRE(I[i][j] == 1.0);
            else REQUIRE(I[i][j] == 0.0);
        }
    }
}
