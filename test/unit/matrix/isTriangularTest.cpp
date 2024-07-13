#include <matrixTest.hpp>

// many of these test can also be done on a large sample size, but three is such a thing as over testing.
TEST_CASE("Does isUpperTriangular() return expected values to corresponding matrices", "[isUpperTriangular]") {
    Matrix<double> upperTriangular;
    Matrix<double> lowerTriangular;
    Matrix<double> sparseUpperMatrix;
    Matrix<double> fullMatrix;
    Matrix<double> identityMatrix;
    upperTriangular.addRow({1, 2, 3})
                   .addRow({0, 4, 5})
                   .addRow({0, 0, 6});
    
    lowerTriangular.addRow({1, 0, 0})
                   .addRow({2, 4, 0})
                   .addRow({3, 5, 6});

    sparseUpperMatrix.addRow({0, 0, 1})
                     .addRow({0, 1, 0})
                     .addRow({0, 0, 0});
    
    fullMatrix.addRow({1, 2, 3})
              .addRow({4, 5, 6})
              .addRow({7, 8, 9});


    REQUIRE(upperTriangular.isUpperTriangular() == true);
    REQUIRE(lowerTriangular.isUpperTriangular() == false);
    REQUIRE(sparseUpperMatrix.isUpperTriangular() == true);
    REQUIRE(fullMatrix.isUpperTriangular() == false);
    
    for (int i = 1; i <= sampleSize; ++i) {
        identityMatrix = Matrix<double>::getIdentity(i);
        REQUIRE(identityMatrix.isUpperTriangular() == true);
    }
}

TEST_CASE("Does isLowerTriangular() return expected values to corresponding matrices", "[isUpperTriangular]") {
    Matrix<double> upperTriangular;
    Matrix<double> lowerTriangular;
    Matrix<double> sparseLowerMatrix;
    Matrix<double> fullMatrix;
    Matrix<double> identityMatrix;
    upperTriangular.addRow({1, 2, 3})
                   .addRow({0, 4, 5})
                   .addRow({0, 0, 6});
    
    lowerTriangular.addRow({1, 0, 0})
                   .addRow({2, 4, 0})
                   .addRow({3, 5, 6});

    sparseLowerMatrix.addRow({0, 0, 0})
                     .addRow({0, 1, 0})
                     .addRow({1, 0, 0});
    
    fullMatrix.addRow({1, 2, 3})
              .addRow({4, 5, 6})
              .addRow({7, 8, 9});


    REQUIRE(upperTriangular.isLowerTriangular() == false);
    REQUIRE(lowerTriangular.isLowerTriangular() == true);
    REQUIRE(sparseLowerMatrix.isLowerTriangular() == true);
    REQUIRE(fullMatrix.isLowerTriangular() == false);
    
    for (int i = 1; i <= sampleSize; ++i) {
        identityMatrix = Matrix<double>::getIdentity(i);
        REQUIRE(identityMatrix.isLowerTriangular() == true);
    }
}

