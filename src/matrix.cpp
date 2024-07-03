#include "matrix.hpp"

template <typename T>
bool Matrix<T>::isMatrixSquare() const {
    return rowSize() > 0 && rowSize() == colSize(); 
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) const {
        if (colSize() != B.rowSize()) {
            std::cerr << "Matrix dimensions do not allow multiplication \n";
            return Matrix<T>();
        }

        Matrix<T> product(rowSize(), B.colSize(), T()); 
        for (int i = 0; i < rowSize(); ++i) {
            for (int j = 0; j < B.colSize(); ++j) {
                T sum = T(); 
                for (int k = 0; k < colSize(); ++k) { 
                    sum += getElement(i, k) * B.getElement(k, j);
                }
                product.assignElement(i, j, sum);
            }
        }
        return product;
    }

template <typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& vec) const {
    const Matrix<T>& A = *this;
    Matrix<T> B;
    B.addCol(vec);
    Matrix<T> productMatrix = (A * B);
    std::vector<T> product = productMatrix.getRow(0);
    return product;
}


template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, const Matrix<T>& B)  {
    Matrix<T> A;
    A.addRow(vec);
    Matrix<T> productMatrix = (A * B);
    std::vector<T> product = productMatrix.getCol(0);
    return product;
}

// For non-const objects
template <typename T>
std::vector<T>& Matrix<T>::operator[](int index) {
    if (index >= rowSize() || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

// For const objects
template <typename T>
const std::vector<T>& Matrix<T>::operator[](int index) const {
    if (index >= rowSize() || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}
template <typename T>
std::vector<T> Matrix<T>::getRow(size_t rowIndex) const {
    std::vector<T> product;
    for (int i = 0; i < rowSize(); ++i) {
        product.push_back(getElement(i, rowIndex));
    }
    return product;
}

template <typename T>
const T& Matrix<T>::getElement(size_t row, size_t col) const {
    return data[row][col];
}

template <typename T>
T& Matrix<T>::getElement(size_t row, size_t col) {
    return data[row][col];
}

template <typename T>
void Matrix<T>::assignElement(size_t row, size_t col, T el) {
    data[row][col] = el;
}

template <typename T>
Matrix<T>& Matrix<T>::addRow(std::vector<T> row) {
    if (row.size() != cols && cols != 0) throw std::runtime_error("Invalid row size");
    if (cols == 0) cols = row.size();
    data.push_back(row);
    rows++;
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::addCol(std::vector<T> col) {
    if (col.size() != rows && rows != 0) throw std::runtime_error("Invalid row size");
    if (rows == 0) rows = col.size();
    for (int i = 0; i < rows; i++) {
        data[i].push_back(col[i]);
    }
    cols++;
    return *this;
}

template <typename T>
void Matrix<T>::print() const {
    for (int i = 0; i < rowSize(); i++) {
        for (int j = 0; j < colSize(); j++) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;