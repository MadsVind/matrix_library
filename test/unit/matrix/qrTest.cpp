#include <matrixTest.hpp>

// Site used for test results:
// https://www.emathhelp.net/en/calculators/linear-algebra/qr-factorization-calculator/?i=%5B%5B1%2C2%2C3%2C4%5D%2C%5B5%2C6%2C7%2C8%5D%2C%5B9%2C10%2C11%2C12%5D%5D

// what happens empty matrix, or is there matrices for which the decomp is not possible?
TEST_CASE("Matrix qr decomposition method", "[qr]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    matrix3.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9})
           .addRow({10, 11, 12});

    SECTION("Normal case with square matrix") {
        Matrix<double>::Qr qr = matrix1.decompQR();

         // Q
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(0), {0.123091490979333,0.904534033733291,0.408248290463863}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(1), {0.492365963917331,0.301511344577764,-0.816496580927726}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(2), {0.861640436855329,-0.301511344577764,0.408248290463863}));
        
        // R
        REQUIRE(checkVecApprox(qr.upper.getRow(0), {8.12403840463596,9.601136296387953,11.078234188139946}));
        REQUIRE(checkVecApprox(qr.upper.getRow(1), {0,0.904534033733291,1.809068067466582}));
        REQUIRE(checkVecApprox(qr.upper.getRow(2), {0, 0, 0, 0}));
        
    }

    SECTION("Normal case with rectangular matrix") {
        Matrix<double>::Qr qr = matrix2.decompQR();

        qr.orthogonal.print();
        std::cout << "\n";
        qr.upper.print();  

        // Q
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(0), {0.096673648904566, 0.907737593658437, 0.408248290463863}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(1), {0.483368244522832, 0.315734815185543, -0.816496580927726}));
        REQUIRE(checkVecApprox(qr.orthogonal.getRow(2), {0.870062840141097, -0.27626796328735, 0.408248290463863}));
        
        // R
        REQUIRE(checkVecApprox(qr.upper.getRow(0), {10.3440804327886, 11.794185166357096, 13.244289899925591, 14.694394633494087}));
        REQUIRE(checkVecApprox(qr.upper.getRow(1), {0, 0.94720444555663, 1.89440889111326, 2.84161333666989}));
        REQUIRE(checkVecApprox(qr.upper.getRow(2), {0, 0, 0, 0}));
        
    }

    SECTION("Normal case with tall matrix") {
        Matrix<double>::Qr qr = matrix3.decompQR();

        qr.orthogonal.print();
        std::cout << "\n";
        qr.upper.print();  
        //REQUIRE(checkVecApprox(refMatrix.getRow(0), {10, 11, 12}));
    }

}