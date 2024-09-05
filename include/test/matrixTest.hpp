#ifndef MATRIX_TEST_HPP
#define MATRIX_TEST_HPP

//#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <matrix.hpp>
#include <limits>

const size_t sampleSize = 300;
const double tolerance = 0.0001;

template <typename T>
bool checkVecApprox(const std::vector<T>& vec, const std::vector<T>& expected) {
    if (vec.size() != expected.size()) {
        return false;
    }
    for (size_t i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i] - expected[i]) > tolerance) {
            return false;
        }
    }
    return true;
}


template <typename T>
bool checkMatrixApprox(const Matrix<T>& mat1, const Matrix<T>& mat2, double tolerance = 1e-4) {
    if (mat1.getRowAmount() != mat2.getRowAmount() || mat1.getColAmount() != mat2.getColAmount()) return false;
    for (size_t i = 0; i < mat1.getRowAmount(); ++i) {
        for (size_t j = 0; j < mat1.getColAmount(); ++j) {
            if (std::abs(mat1[i][j] - mat2[i][j]) > tolerance) return false;
        }
    }
    return true;
}

template <typename T>
Matrix<T> createIncrementingSquareMatrix(size_t size) {
    Matrix<T> matrix;
    for (int i = 0; i < size; ++i) {
        std::vector<T> row;
        for (int j = 0; j < size; ++j) {
            row.push_back(1 + i * size + j);
        }
        matrix.addRow(row);
    }
    return matrix;
}


#endif