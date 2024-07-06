#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

template <typename T>
class Matrix {

private:
    /* data */
    std::vector<std::vector<T>> data;
    int rowAmount;
    int colAmount;
    
public:
    Matrix() : rowAmount(0), colAmount(0) {}
    Matrix(int rowAmount, int colAmount, T value) : rowAmount(rowAmount), colAmount(colAmount) {
        data.resize(rowAmount);
        for (int i = 0; i < rowAmount; i++) {
            data[i].resize(colAmount, value);
        }
    }

    Matrix(int rowAmount, int colAmount) : rowAmount(rowAmount), colAmount(colAmount) {
        data.resize(rowAmount);
        for (int i = 0; i < rowAmount; i++) {
            data[i].resize(colAmount);
        }
    }

    Matrix(const Matrix<T>& x) : rowAmount(x.rowAmount), colAmount(x.colAmount), data(x.data) {}

    class Plu {
       public: 
        Plu(Matrix<T> permutation, Matrix<T> lower, Matrix<T> upper, size_t permutationAmount) 
            : permutation(permutation), lower(lower), upper(upper), permutationAmount(permutationAmount) {}
        const Matrix<T> permutation;
        const Matrix<T> lower;
        const Matrix<T> upper;
        const size_t permutationAmount;
    };

    static Matrix<T> getIdentity(size_t size);

    size_t getRowAmount() const { return rowAmount; }

    size_t getColAmount() const { return colAmount; }

    bool operator==(const Matrix<T>& B) const;

    Matrix<T> operator-(const Matrix<T>& B) const;
    Matrix<T> operator+(const Matrix<T>& B) const;

    Matrix<T> pow(const T& exponent) const;
    
    Matrix<T> operator*(const Matrix<T>& B) const;
    Matrix<T> operator*(const T& scalar) const;
    std::vector<T> operator*(const std::vector<T>& vec) const;

    // Unary
    T determinant() const;
    Matrix<T> inverse() const;
    Matrix<T> transpose() const;
    Matrix<T> operator-() const;

    Matrix<T>::Plu decompPLU() const;

    std::vector<T>& operator[](int index);

    // Const overload of [] operator for const objects
    const std::vector<T>& operator[](int index) const;

    std::vector<T> getRow(size_t rowIndex) const;
    std::vector<T> getCol(size_t colIndex) const;

    const T& getElement(size_t row, size_t col) const;

    T& getElement(size_t row, size_t col);

    void assignElement(size_t row, size_t col, T el);
    
    Matrix<T>& addRow(std::vector<T> row);
    Matrix<T>& addCol(std::vector<T> col);

    Matrix<T>& setRow(size_t rowIndex, std::vector<T> row);
    Matrix<T>& setCol(size_t colIndex, std::vector<T> col);

    Matrix<T>& swapRow(size_t rowA, size_t rowB);
    Matrix<T>& swapCol(size_t colA, size_t colB);

    void print() const;
    
    bool isMatrixSquare() const;

    bool areMatricesSameSize(const Matrix<T>& other) const;

};

template <typename T>
std::vector<T>& operator*(const std::vector<T>& vec, const Matrix<T>& B); // !!! DOES NOT WORK

template <typename T>
Matrix<T>  operator*(const T& scalar, const Matrix<T>& B); // !!! DOES NOT WORK

template <typename T>
Matrix<T> operator%(const T& scalar, const Matrix<T>& B); // !!! DOES NOT WORK
#endif
