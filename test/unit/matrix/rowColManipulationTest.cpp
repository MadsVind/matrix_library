#include <matrixTest.hpp>

TEST_CASE("Matrix row operations", "[unit_test]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    SECTION("getRow") {
        REQUIRE(matrix1.getRow(0) == std::vector<double>({1, 2, 3}));
        REQUIRE_THROWS_AS(matrix1.getRow(3), std::out_of_range);
    }

    SECTION("addRow") {
        matrix1.addRow({10, 11, 12});
        REQUIRE(matrix1.getRow(3) == std::vector<double>({10, 11, 12}));
        REQUIRE_THROWS_AS(matrix1.addRow({1, 2}), std::runtime_error);
    }

    SECTION("setRow") {
        matrix1.setRow(0, {10, 11, 12});
        REQUIRE(matrix1.getRow(0) == std::vector<double>({10, 11, 12}));
        REQUIRE_THROWS_AS(matrix1.setRow(3, {1, 2, 3}), std::runtime_error);
        REQUIRE_THROWS_AS(matrix1.setRow(0, {1, 2}), std::runtime_error);
    }

    SECTION("swapRow") {
        matrix1.swapRow(0, 1);
        REQUIRE(matrix1.getRow(0) == std::vector<double>({4, 5, 6}));
        REQUIRE(matrix1.getRow(1) == std::vector<double>({1, 2, 3}));
        REQUIRE_THROWS_AS(matrix1.swapRow(3, 0), std::runtime_error);
        REQUIRE_THROWS_AS(matrix1.swapRow(0, 3), std::runtime_error);
    }
}

TEST_CASE("Matrix col operations", "[col operations]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});
    
    SECTION("getCol") {
        REQUIRE(matrix1.getCol(0) == std::vector<double>({1, 4, 7}));
        REQUIRE_THROWS_AS(matrix1.getCol(3), std::out_of_range);
    }

    SECTION("addCol") {
        matrix1.addCol({10, 11, 12});
        REQUIRE(matrix1.getCol(3) == std::vector<double>({10, 11, 12}));
        REQUIRE_THROWS_AS(matrix1.addCol({1, 2}), std::runtime_error);
    }

    SECTION("setCol") {
        matrix1.setCol(0, {10, 11, 12});
        REQUIRE(matrix1.getCol(0) == std::vector<double>({10, 11, 12}));
        REQUIRE_THROWS_AS(matrix1.setCol(3, {1, 2, 3}), std::runtime_error);
        REQUIRE_THROWS_AS(matrix1.setCol(0, {1, 2}), std::runtime_error);
    }

    SECTION("swapCol") {
        matrix1.swapCol(0, 1);
        REQUIRE(matrix1.getCol(0) == std::vector<double>({2, 5, 8}));
        REQUIRE(matrix1.getCol(1) == std::vector<double>({1, 4, 7}));
        REQUIRE_THROWS_AS(matrix1.swapCol(3, 0), std::runtime_error);
        REQUIRE_THROWS_AS(matrix1.swapCol(0, 3), std::runtime_error);
    }
}

TEST_CASE("Matrix element operations", "[element operations]") {
    Matrix<double> matrix1;
    Matrix<double> matrix2;

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    SECTION("operator[]") {
        REQUIRE(matrix1[0] == std::vector<double>({1, 2, 3}));
        REQUIRE_THROWS_AS(matrix1[-1], std::out_of_range);
        REQUIRE_THROWS_AS(matrix1[3], std::out_of_range);
    }
    
    SECTION("getElement") {
        REQUIRE(matrix1.getElement(0, 0) == 1);
        REQUIRE_THROWS_AS(matrix1.getElement(3, 0), std::out_of_range);
        REQUIRE_THROWS_AS(matrix1.getElement(0, 3), std::out_of_range);
    }

    SECTION("assignElement") {
        matrix1.assignElement(0, 0, 10);
        REQUIRE(matrix1.getElement(0, 0) == 10);
        REQUIRE_THROWS_AS(matrix1.assignElement(3, 0, 10), std::out_of_range);
        REQUIRE_THROWS_AS(matrix1.assignElement(0, 3, 10), std::out_of_range);
    }
}