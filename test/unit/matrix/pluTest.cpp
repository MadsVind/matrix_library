#include <matrixTest.hpp>

// Test cases for decompPLU method
TEST_CASE("Matrix plu method", "[unit_test]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;

    matrix1.addRow({2, 0, 2, 0.6})
           .addRow({3, 3, 4, -2})
           .addRow({5, 5, 4, 2})
           .addRow({-1, -2, 3.4, -1});

    matrix2.addRow({1, 0, 0})
           .addRow({0, 1, 0})
           .addRow({0, 0, 1});

    matrix3.addRow({0, 2, 2})
           .addRow({1, 2, 3})
           .addRow({4, 4, 4});

    SECTION("Simple matrix") {
        Matrix<double>::Plu lu = matrix1.plu();
        // Expected permutation matrix
        Matrix<double> expectedPerm;
        expectedPerm.addRow({0, 0, 1, 0})
                    .addRow({1, 0, 0, 0})
                    .addRow({0, 0, 0, 1})
                    .addRow({0, 1, 0, 0});
        // Expected lower matrix
        Matrix<double> expectedLower;
        expectedLower.addRow({1,    0,   0,   0})
                     .addRow({0.4,  1,   0,   0})
                     .addRow({-0.2, 0.5, 1,   0})
                     .addRow({0.6,  0,   0.4, 1});

        // Expected upper matrix
        Matrix<double> expectedUpper;
        expectedUpper.addRow({5,  5, 4,     2})
                     .addRow({0, -2, 0.4,  -0.2})
                     .addRow({0,  0, 4,    -0.5})
                     .addRow({0,  0, 0,    -3});

        // Compare permutation, lower, and upper matrices
        REQUIRE(checkMatrixApprox(lu.permutation, expectedPerm));
        REQUIRE(checkMatrixApprox(lu.lower, expectedLower));
        REQUIRE(checkMatrixApprox(lu.upper, expectedUpper));
    }

    SECTION("Identity matrix") {
        Matrix<double>::Plu lu = matrix2.plu();

        // Expected permutation, lower, and upper matrices are identity matrices
        Matrix<double> expectedIdentity;
        expectedIdentity.addRow({1, 0, 0})
                        .addRow({0, 1, 0})
                        .addRow({0, 0, 1});

        // Compare permutation, lower, and upper matrices
        REQUIRE(checkMatrixApprox(lu.permutation, expectedIdentity));
        REQUIRE(checkMatrixApprox(lu.lower, expectedIdentity));
        REQUIRE(checkMatrixApprox(lu.upper, expectedIdentity));
    }

    SECTION("Matrix requiring row swaps") {
        Matrix<double>::Plu lu = matrix3.plu();

        // Expected permutation matrix
        Matrix<double> expectedPerm;
        expectedPerm.addRow({0, 0, 1})
                    .addRow({1, 0, 0})
                    .addRow({0, 1, 0});

        // Expected lower matrix
        Matrix<double> expectedLower;
        expectedLower.addRow({1, 0, 0})
                     .addRow({0, 1, 0})
                     .addRow({0.25, 0.5, 1});

        // Expected upper matrix
        Matrix<double> expectedUpper;
        expectedUpper.addRow({4, 4, 4})
                     .addRow({0, 2, 2})
                     .addRow({0, 0, 1});

        // Compare permutation, lower, and upper matrices
        REQUIRE(checkMatrixApprox(lu.permutation, expectedPerm));
        REQUIRE(checkMatrixApprox(lu.lower, expectedLower));
        REQUIRE(checkMatrixApprox(lu.upper, expectedUpper));
    }
}