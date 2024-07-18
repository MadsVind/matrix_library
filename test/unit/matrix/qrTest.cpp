#include <matrixTest.hpp>

// Site used for test results:
// https://www.emathhelp.net/en/calculators/linear-algebra/qr-factorization-calculator/?i=%5B%5B1%2C2%2C3%2C4%5D%2C%5B5%2C6%2C7%2C8%5D%2C%5B9%2C10%2C11%2C12%5D%5D

// what happens empty matrix, or is there matrices for which the decomp is not possible?
TEST_CASE("Matrix qr decomposition method", "[qr]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;

    matrix1.addRow({1, 1, 0})
           .addRow({1, 0, 1})
           .addRow({0, 1, 1});

    matrix2.addRow({1, 1, 0, 1})
           .addRow({1, 0, 1, 1})
           .addRow({0, 1, 1, 1});

    matrix3.addRow({1, 1, 0})
           .addRow({1, 0, 1})
           .addRow({0, 1, 1})
           .addRow({1, 0, 1});

    SECTION("Normal case with square matrix") {
        Matrix<double>::Qr qr = matrix1.decompQR();

         // Q
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(0), {0.707107, 0.408248, -0.57735}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(1), {0.707107, -0.408248, 0.57735}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(2), {0, 0.816497, 0.57735 }));
        
        // R
        REQUIRE(checkVecApprox(qr.upper.getRow(0), {1.41421, 0.707107, 0.707107 }));
        REQUIRE(checkVecApprox(qr.upper.getRow(1), {0, 1.22474, 0.408248}));
        REQUIRE(checkVecApprox(qr.upper.getRow(2), {0, 0, 1.1547}));
        
    }

    SECTION("Normal case with rectangular matrix") {
        Matrix<double>::Qr qr = matrix2.decompQR();

        // Q
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(0), {0.707107, 0.408248, -0.57735}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(1), {0.707107, -0.408248, 0.57735}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(2), {0, 0.816497, 0.57735}));
        
        // R
        REQUIRE(checkVecApprox(qr.upper.getRow(0), {1.41421, 0.707107, 0.707107, 1.41421}));
        REQUIRE(checkVecApprox(qr.upper.getRow(1), {0, 1.22474, 0.408248, 0.816497}));
        REQUIRE(checkVecApprox(qr.upper.getRow(2), {0, 0, 1.1547, 0.57735}));
        
    }

    SECTION("Normal case with tall matrix") {
        Matrix<double>::Qr qr = matrix3.decompQR();

        // Q
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(0), {0.57735, 0.516398, -0.632456, 0}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(1), {0.57735, -0.258199, 0.316228, 0}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(2), {0, 0.774597, 0.632456, 0 }));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(3), {0.57735, -0.258199, 0.316228, 0}));
        
        // R
        REQUIRE(checkVecApprox(qr.upper.getRow(0), {1.73205, 0.57735, 1.1547}));
        REQUIRE(checkVecApprox(qr.upper.getRow(1), {0, 1.29099, 0.258199}));
        REQUIRE(checkVecApprox(qr.upper.getRow(2), {0, 0, 1.26491}));
        REQUIRE(checkVecApprox(qr.upper.getRow(3), {0, 0, 0}));
    }

}