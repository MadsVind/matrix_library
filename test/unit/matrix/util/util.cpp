#include <matrixTest.hpp>

bool checkVecApprox(const std::vector<double>& vec, const std::vector<double>& expected) {
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