#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "matrix.hpp"

template <typename T>
Matrix<T> parseMatrix(const std::string& str) {
    std::istringstream iss(str);
    std::vector<std::vector<T>> data;
    std::string line;
    while (std::getline(iss, line, ';')) {
        std::istringstream lineStream(line);
        std::vector<T> row;
        T value;
        while (lineStream >> value) {
            row.push_back(value);
        }
        data.push_back(row);
    }
    int rows = data.size();
    int cols = data[0].size();
    Matrix<T> matrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix.assignElement(i, j, data[i][j]);
        }
    }
    return matrix;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <method> <matrix1> [<matrix2>]" << std::endl;
        return 1;
    }

    std::string method = argv[1];
    Matrix<double> matrix1 = parseMatrix<double>(argv[2]);
    Matrix<double> matrix2;
    if (argc == 4) {
        matrix2 = parseMatrix<double>(argv[3]);
    }

    if (method == "add" && argc == 4) {
        Matrix<double> result = matrix1 + matrix2;
        result.print();
    } else if (method == "subtract" && argc == 4) {
        Matrix<double> result = matrix1 - matrix2;
        result.print();
    } else if (method == "dot" && argc == 4) {
        Matrix<double> result = matrix1 * matrix2;
        result.print();
    } else if (method == "pow" && argc == 4) {
        Matrix<double> result = matrix1.pow(matrix2[0][0]);
        result.print();
    } else if (method == "determinant") {
        double result = matrix1.determinant();
        std::cout << "Determinant: " << result << std::endl;
    } else if (method == "inverse") {
        Matrix<double> result = matrix1.inverse();
        result.print();
    } else if (method == "transpose") {
        Matrix<double> result = matrix1.transpose();
        result.print();
    } else if (method == "ref") {
        Matrix<double> result = matrix1.ref();
        result.print();
    } else if (method == "rref") {
        Matrix<double> result = matrix1.rref();
        result.print();
    } else if (method == "eigen") {
        Matrix<double>::Eigen result = matrix1.eigen();
        result.vectorVec.print();
        std::cout << "\n";
        for (const auto& val : result.valueVec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    } else if (method == "qr") {
        Matrix<double>::Qr result = matrix1.qr();
        result.orthogonal.print();
        std::cout << "\n";
        result.upper.print();
    } else if (method == "plu") {
        Matrix<double>::Plu result = matrix1.plu();
        result.permutation.print();
        std::cout << "\n";
        result.lower.print();
        std::cout << "\n";
        result.upper.print();
    } else {
        std::cerr << "Unknown method or incorrect number of arguments" << std::endl;
        return 1;
    }

    return 0;
}