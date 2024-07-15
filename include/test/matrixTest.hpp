#ifndef MATRIX_TEST_HPP
#define MATRIX_TEST_HPP

#include <catch2/catch_test_macros.hpp>
#include <matrix.hpp>
#include <limits>

const size_t sampleSize = 1000;
const double tolerance = 0.0001;

bool checkVecApprox(const std::vector<double>& vec, const std::vector<double>& expected);

#endif