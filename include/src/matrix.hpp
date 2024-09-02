#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>
#include <thread>
#include <mutex>

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

    Matrix(const Matrix<T>& x) : rowAmount(x.rowAmount), colAmount(x.colAmount), data(x.data.size()) {
        for (size_t i = 0; i < x.data.size(); ++i) {
            data[i] = x.data[i];
        }
    }

    class Plu {
       public: 
        Plu(const Matrix<T>& permutation, const Matrix<T>& lower, const Matrix<T>& upper, const size_t& permutationAmount) 
            : permutation(permutation), lower(lower), upper(upper), permutationAmount(permutationAmount) {}
        const Matrix<T> permutation;
        const Matrix<T> lower;
        const Matrix<T> upper;
        const size_t permutationAmount;
    };

    class Eigen {
       public: 
        Eigen(const Matrix<T>& vectorVec, const std::vector<T>& valueVec) 
            : vectorVec(vectorVec), valueVec(valueVec) {}
        const Matrix<T> vectorVec;
        const std::vector<T> valueVec;
    };

    class Qr {
     public: 
        Qr(const Matrix<T>& orthogonal, const Matrix<T>& upper) 
            : orthogonal(orthogonal), upper(upper) {}
        const Matrix<T> orthogonal;
        const Matrix<T> upper;
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

    Matrix<T> ref(double tolerance = std::numeric_limits<T>::min()) const;
    Matrix<T> rref(double tolerance = std::numeric_limits<T>::min()) const;
    Matrix<T>::Eigen eigen(double tolerance = 0.00000000000000000001) const; //! magick number for tolerance to work apperently. !!! Test this
    Matrix<T>::Qr qr() const;
    Matrix<T>::Qr givensRotations() const;
    Matrix<T>::Qr gramSmith() const;
    Matrix<T>::Plu plu() const;

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

    void print(double tolerance = 0) const;
    
    bool isSquare() const;

    bool isSameSize(const Matrix<T>& other) const;

    bool isUpperTriangular(double tolerance = std::numeric_limits<T>::min()) const;

    bool isLowerTriangular(double tolerance = std::numeric_limits<T>::min()) const;

};

template <typename T>
std::vector<T>& operator*(const std::vector<T>& vec, const Matrix<T>& B); // !!! DOES NOT WORK

template <typename T>
Matrix<T>  operator*(const T& scalar, const Matrix<T>& B); // !!! DOES NOT WORK

template <typename T>
Matrix<T> operator%(const T& scalar, const Matrix<T>& B); // !!! DOES NOT WORK
#endif
