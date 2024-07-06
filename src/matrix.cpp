#include "matrix.hpp"


template <typename T>
Matrix<T> Matrix<T>::getIdentity(size_t size) {
    Matrix<T> identity(size, size, static_cast<T>(0));
    T one = static_cast<T>(1); // maybe try catch this
    for (int i = 0; i < size; ++i) {
        identity[i][i] = one;
    }
    return identity;
}

template <typename T>
bool Matrix<T>::isMatrixSquare() const {
    return getRowAmount() > 0 && getRowAmount() == getColAmount(); 
}

template <typename T>
bool Matrix<T>::areMatricesSameSize(const Matrix<T>& B) const {
    return getRowAmount() == B.getRowAmount() && getColAmount() == B.getColAmount();
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& B) const {
    if (!areMatricesSameSize(B)) return false;
    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();

    Matrix<T> product(rowAmount, colAmount);
    for (int i = 0; i < rowAmount ; ++i) {
        for (int j = 0; j <colAmount; ++j) {
            if (data[i][j] != B[i][j]) return false;
        }
    }
    return true;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& B) const {
    if (!areMatricesSameSize(B)) {
        std::cerr << "Matrices was not of matching sizes and can therefore no be added together";
        return Matrix<T>();
    }
    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    Matrix<T> product(rowAmount , colAmount);

    for (int i = 0; i < rowAmount ; ++i) {
        for (int j = 0; j < colAmount; ++j) {
            product[i][j] = data[i][j] + B[i][j];
        }
    }
    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& B) const {
    if (!areMatricesSameSize(B)) {
        std::cerr << "Matrices was not of matching sizes and can therefore no be added together";
        return Matrix<T>();
    }
    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    Matrix<T> product(rowAmount , colAmount);

    for (int i = 0; i < rowAmount; ++i) {
        for (int j = 0; j < colAmount; ++j) {
            product[i][j] = data[i][j] - B[i][j];
        }
    }
    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const {
    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    Matrix<T> product(rowAmount, colAmount);

    for (int i = 0; i < rowAmount; ++i) {
        for (int j = 0; j < colAmount; ++j) {
            product[i][j] = - data[i][j];
        }
    }
    return product;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    Matrix<T> product;

    for (int i = 0; i < rowAmount; ++i) {
        std::vector<T> row = getRow(i);
        product.addCol(row);
    }
    return product;
}



template <typename T>
typename Matrix<T>::Plu Matrix<T>::decompPLU() const {
    Matrix<T> perm = getIdentity(rowAmount);
    Matrix<T> lower = Matrix<T>(perm);
    Matrix<T> upper = Matrix<T>(*this);

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    
    size_t permutationAmount = 0;
    int pivotCol = 0;
    for (int i = 0; i < rowAmount; ++i) { // Find pivot in every row
        int largestPivotRow = -1;
        
        for (int j = pivotCol; j < colAmount; ++j) { // Go to next column if column has no pivot
            for (int k = i; k < rowAmount; ++k) { 

                T currentElement = upper[k][j];
                if ((largestPivotRow == -1 && currentElement != 0) || 
                    (largestPivotRow != -1 && currentElement > std::abs(upper[largestPivotRow][j]))) {
                        largestPivotRow = k;
                }
            }
            pivotCol = j;
            if (largestPivotRow != -1) break; // no pivot only zeroes if this is the case
        }
        if (pivotCol > rowAmount) break;
        
        if (largestPivotRow != i) {
            ++permutationAmount;
            upper.swapRow(i, largestPivotRow);
            perm.swapRow(i, largestPivotRow);
            lower.swapRow(i, largestPivotRow);
            lower.swapCol(i, largestPivotRow);  // this is a thing which is done :D i don't get why.
        }

        T zero = static_cast<T>(0);
        for (int k = i + 1; k < rowAmount; ++k) {
            if (upper[k][pivotCol] == 0) continue;
            T scalar = upper[k][pivotCol] / upper[i][pivotCol]; 
            lower[k][pivotCol] = scalar;  // sometimes this should be minus other time not?!?!?
            for (int l = 0; l < colAmount; ++l) {
                upper[k][l] -= scalar * upper[i][l];
            }
        }
    }
    return Matrix<T>::Plu(perm, lower, upper, permutationAmount);
}



template <typename T>
T Matrix<T>::determinant() const {
    Matrix<T>::Plu plu = decompPLU();
    Matrix<T> upper = plu.upper;
    size_t size = upper.getRowAmount();
    T determinant = upper[0][0];
    for (int i = 1; i < size; ++i) {
        determinant *= upper[i][i];
    }
    if (plu.permutationAmount % 2 != 0) determinant *= -1;
    return determinant;
}

// get 1 then clear column (ex. get 1 in in 1st row 1st col, now clear rest of 1st col to 0, repeat for 2nd row 2nd col ect.)
template <typename T>
Matrix<T> Matrix<T>::inverse() const {
    if (!isMatrixSquare()) {
        std::cerr << "Can not get inverse of non square matrix \n";
        return Matrix<T>();
    }

    if (determinant() == 0) {
        std::cerr << "Can not get inverse of matrix with determinant of 0 \n";
        return Matrix<T>();
    }

    size_t size = getRowAmount();
    Matrix<T> iden = getIdentity(size);
    Matrix<T> org(*this);


    for (int i = 0; i < size; ++i) {
        T scalar = org[i][i];
        for (int j = 0; j < size; ++j) { 
            org[i][j] = org[i][j] / scalar;
            iden[i][j] = iden[i][j] / scalar;
        }

        for (int j = 0; j < size; ++j) {
            if (i == j) continue;
            T scalar = org[j][i]; 
            for (int k = 0; k < size; ++k) {
                org[j][k] -= scalar * org[i][k];
                iden[j][k] -= scalar * iden[i][k];
            }
        }
    }
    return iden;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) const {
    if (getColAmount() != B.getRowAmount()) {
        std::cerr << "Matrix dimensions do not allow multiplication \n";
        return Matrix<T>();
    }

    Matrix<T> product(getRowAmount(), B.getColAmount(), T()); 
    for (int i = 0; i < getRowAmount(); ++i) {
        for (int j = 0; j < B.getColAmount(); ++j) {
            T sum = T(); 
            for (int k = 0; k < getColAmount(); ++k) { 
                sum += data[k][j] * B[i][k];
            }
            product[i][j] = sum;
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
Matrix<T> Matrix<T>::operator*(const T& scalar) const  {
    size_t colAmount = getColAmount();
    size_t rowAmount = getRowAmount();

    Matrix<T> product(rowAmount, colAmount);
    for (size_t i = 0; i < rowAmount; ++i) {
        for (size_t j = 0; j < colAmount; ++j) {
            product[i][j] = scalar * getElement(i, j); 
        }
    }
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

template <typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& B)  {
    return B * scalar;
}

//Matrix<T> pow(const T& exponent) const;


// For non-const objects
template <typename T>
std::vector<T>& Matrix<T>::operator[](int index) {
    if (index >= getRowAmount() || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

// For const objects
template <typename T>
const std::vector<T>& Matrix<T>::operator[](int index) const {
    if (index >= getRowAmount() || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

template <typename T>
std::vector<T> Matrix<T>::getCol(size_t colIndex) const {
    std::vector<T> product;
    for (int i = 0; i < getRowAmount(); ++i) {
        product.push_back(data[i][colIndex]);
    } 
    return product;
}

template <typename T>
std::vector<T> Matrix<T>::getRow(size_t rowIndex) const {
    return data[rowIndex];
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
    if (row.size() != colAmount && colAmount != 0) throw std::runtime_error("Invalid row size");
    if (colAmount == 0) colAmount = row.size();
    data.push_back(row);
    rowAmount++;
    return *this;
}


template <typename T>
Matrix<T>& Matrix<T>::addCol(std::vector<T> col) {
    if (col.size() != rowAmount && rowAmount != 0) throw std::runtime_error("Invalid row size");
    bool isEmpty = rowAmount == 0;
    if (isEmpty) rowAmount = col.size();
    for (int i = 0; i < rowAmount; i++) {
        if (isEmpty) data.push_back(std::vector<T>());
        data[i].push_back(col[i]);
    }
    colAmount++;
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::setRow(size_t rowIndex, std::vector<T> row) {
    if (rowIndex >= rowAmount) throw std::runtime_error("Tried to set non existing row");
    if (row.size() != colAmount) throw std::runtime_error("Invalid row size");
    for (int i = 0; i < colAmount; ++i) {
        data[rowIndex][i] = row[i];
    }

    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::setCol(size_t colIndex, std::vector<T> col) {
    if (colIndex >= colAmount) throw std::runtime_error("Tried to set non existing col");
    if (col.size() != rowAmount) throw std::runtime_error("Invalid col size");
    for (int i = 0; i < rowAmount; ++i) {
        data[i][colIndex] = col[i];
    }
    return *this;
}

// if you swap the same row with the same row you deserve bad performance
template <typename T>
Matrix<T>& Matrix<T>::swapRow(size_t rowA, size_t rowB) { 
    if (rowA > rowAmount || rowB > rowAmount) throw std::runtime_error("Tried to swap non existing row");
    std::vector<T> tempRow = data[rowA];
    (*this).setRow(rowA, data[rowB]);
    (*this).setRow(rowB, tempRow);
    return *this;
}

template <typename T>
 Matrix<T>& Matrix<T>::swapCol(size_t colA, size_t colB) {
     if (colA > colAmount || colB > colAmount) throw std::runtime_error("Tried to swap non existing col");
    std::vector<T> tempCol = (*this).getCol(colA);
    (*this).setCol(colA, (*this).getCol(colB));
    (*this).setCol(colB, tempCol);
    return *this;
}

template <typename T>
void Matrix<T>::print() const {
    for (int i = 0; i < getRowAmount(); i++) {
        for (int j = 0; j < getColAmount(); j++) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;