#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>

template <typename T>
class Matrix {

private:
    /* data */
    std::vector<std::vector<T>> data;
    int rows;
    int cols;
    
public:
    Matrix() : rows(0), cols(0) {}
    Matrix(int rows, int cols, T value) : rows(rows), cols(cols) {
        data.resize(rows);
        for (int i = 0; i < rows; i++)
        {
            data[i].resize(cols, value);
        }
    }

    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data.resize(rows);
        for (int i = 0; i < rows; i++)
        {
            data[i].resize(cols);
        }
    }

    size_t rowSize() const { return rows; }

    size_t colSize() const { return cols; }

    bool operator==(const Matrix<T>& B) const;

    Matrix<T> operator*(const Matrix<T>& B) const;

    std::vector<T> operator*(const std::vector<T>& vec) const;
    
    Matrix<T> operator*(const T& scalar) const;

    std::vector<T>& operator[](int index);

    // Const overload of [] operator for const objects
    const std::vector<T>& operator[](int index) const;

    std::vector<T> getRow(size_t rowIndex) const;

    const T& getElement(size_t row, size_t col) const;

    T& getElement(size_t row, size_t col);

    void assignElement(size_t row, size_t col, T el);
    
    Matrix<T>& addRow(std::vector<T> row);

    Matrix<T>& addCol(std::vector<T> col);

    void print() const;
    
    bool isMatrixSquare() const;
};

template <typename T>
std::vector<T>& operator*(const std::vector<T>& vec, const Matrix<T>& B); // !!! DOES NOT WORK

template <typename T>
Matrix<T>  operator*(const T& scalar, const Matrix<T>& B); // !!! DOES NOT WORK

#endif