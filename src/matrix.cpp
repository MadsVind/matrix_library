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
bool Matrix<T>::isSquare() const {
    return getRowAmount() > 0 && getRowAmount() == getColAmount(); 
}

template <typename T>
bool Matrix<T>::isSameSize(const Matrix<T>& B) const {
    return getRowAmount() == B.getRowAmount() && getColAmount() == B.getColAmount();
}

template <typename T>
bool Matrix<T>::isUpperTriangular(double tolerance) const {
    if (!isSquare()) return false;
    size_t size = getRowAmount();
    for (size_t i = 1; i < size; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::abs(data[i][j]) > tolerance) return false;
        }
    }
    return true;
}

template <typename T>
bool Matrix<T>::isLowerTriangular(double tolerance) const {
    if (!isSquare()) return false;
    size_t size = getRowAmount();
    for (size_t i = 1; i < size; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::abs(data[j][i]) > tolerance) return false;
        }
    }
    return true;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& B) const {
    if (!isSameSize(B)) return false;
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
    if (!isSameSize(B)) throw std::invalid_argument("Matrices was not of matching sizes and can therefore no be added together");
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
    if (!isSameSize(B)) throw std::invalid_argument("Matrices was not of matching sizes and can therefore no be added together");
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
Matrix<T> Matrix<T>::ref(double tolerance) const {
    Matrix<T> ref(*this);

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    
    T zero = static_cast<T>(0);

    int pivot = 0;
    for (int col = 0; col < colAmount; ++col) {
        int largestPivotRowIndex = -1;
        
        // Find largest pivot in column
        for (int tempRow = pivot; tempRow < rowAmount; ++tempRow) {          
            T currentElement = std::abs(ref[tempRow][col]);
            if (largestPivotRowIndex == -1 && currentElement > tolerance ||
                largestPivotRowIndex != -1 && currentElement > std::abs(ref[largestPivotRowIndex][col])) {
                largestPivotRowIndex = tempRow; // Should work since 
            }
        }

        // If no pivot in column then continue
        if (largestPivotRowIndex == -1) continue;
        // Swap rows if standard pivot wasn't the abs largest
        if (largestPivotRowIndex != pivot) ref.swapRow(pivot, largestPivotRowIndex);
        
        // Eliminate values lower in pivot column
        for (int row = pivot + 1; row < rowAmount; ++row) {
            T nonPivotInCol = ref[row][col];
            if (std::abs(nonPivotInCol) < tolerance) continue;

            T scalar = nonPivotInCol / ref[pivot][col]; 

            // Minus pivot row from lower row 
            for (int i = col; i < colAmount; ++i) {
                ref[row][i] -= scalar * ref[pivot][i];
            }
        }
        ++pivot;
    }
    return ref;
}


template <typename T>
Matrix<T> Matrix<T>::rref(double tolerance) const { 
    Matrix<T> rref(*this);

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    
    T zero = static_cast<T>(0);

    int pivot = 0;
    for (int col = 0; col < colAmount; ++col) {
        int largestPivotRowIndex = -1;
        
        for (int tempRow = pivot; tempRow < rowAmount; ++tempRow) {       
            T currentElement = std::abs(rref[tempRow][col]);
            if (largestPivotRowIndex == -1 && currentElement > tolerance ||
                largestPivotRowIndex != -1 && currentElement > std::abs(rref[largestPivotRowIndex][col])) {
                largestPivotRowIndex = tempRow;
            }
        }

        if (largestPivotRowIndex == -1) continue;
        if (largestPivotRowIndex != pivot) rref.swapRow(pivot, largestPivotRowIndex);

        T scalar = rref[pivot][col];

        for (int i = col; i < colAmount; ++i) { 
            rref[pivot][i] = rref[pivot][i] / scalar;
        }
    
        for (int row = 0; row < rowAmount; ++row) {
            T nonPivotInCol = rref[row][col];
            if (pivot == row || nonPivotInCol == 0) continue;

            T scalar = nonPivotInCol / rref[pivot][col]; 
            for (int i = col; i < colAmount; ++i) {
                rref[row][i] -= scalar * rref[pivot][i];
            }
        }
        ++pivot;
    }
    return rref;
}

template <typename T>
std::vector<T> backSubstitution(Matrix<T> A, double tolerance) { // Doesn't work if pivot is not on diagonal
    size_t rowAmount = A.getRowAmount();
    size_t colAmount = A.getColAmount();

    size_t rowAmountIndex = (rowAmount - 1);

    std::vector<T> res(rowAmount);
    T one = static_cast<T>(1);
    T zero = static_cast<T>(0);
    for (int row = rowAmountIndex; 0 <= row; --row) {
        bool isEmpty = true;
        int pivot = -1;
        for (int col = row; col < colAmount; ++col) { 
            if (std::abs(A[row][col]) > tolerance) {
                isEmpty = false;
                pivot = col;
                break;
            }
        }

        if (isEmpty) {
            res[row] = one; 
            continue;
        }

        T varSum = zero;
        T pivotVal = A[row][pivot];

        int pivotOffset = row - pivot;
        for (int col = pivot + 1; col < colAmount; ++col) {
           varSum -= (A[row][col] * res[col + pivotOffset]) / pivotVal; 
        }
        res[row] = varSum;
    }
    std::cout << "\n";
    return res;
}

template <typename T>
typename Matrix<T>::Qr Matrix<T>::decompQR() const {
    if (determinant() == 0) {
        std::cerr << "Can not use Gram-Schmidt method on matrix since it is linuarly dependent\n";
        return Matrix<T>::Qr(Matrix<T>(), Matrix<T>()); 
    }
    T zero = static_cast<T>(0);
    Matrix<T> orthogonal = Matrix<T>(rowAmount, rowAmount);
    Matrix<T> upper = Matrix<T>(rowAmount, colAmount, zero);

    // For every column of A
    for (int col = 0; col < colAmount; ++col) {
        std::vector<T> a = getCol(col);
        std::vector<T> scalars;

        // Calculate the scalars of the new orthogonal vector a (alpha)
        // These scalars are also the non digonal entries in upper
        for (int qIdx = 0; qIdx < col; ++qIdx) {
            T scalar = zero;
            std::vector<T> q = orthogonal.getCol(qIdx);

            // Vector product of columns of org matrix A and Q, which will become a scalar 
            for (int row = 0; row < a.size(); ++row) {
                scalar += a[row] * q[row];
            }
            scalars.push_back(scalar);
            upper[qIdx][col] = scalar; 
        }

        // Calculate the new orthogonal vector a (alpha)
        std::vector<T> alpha(a.size(), zero);
        for (int row = 0; row < a.size(); ++row) {
            T el = a[row];
            for (int qIdx = 0; qIdx < scalars.size(); ++qIdx) {
                std::vector<T> q = orthogonal.getCol(qIdx);
                el -= scalars[qIdx] * q[row];
            }
            alpha[row] = el;
        }
        
        //  Calculate the length of new orthogonal vector a (alpha)
        T alphaLength = zero;
        for (int row = 0; row < alpha.size(); ++row) {
            alphaLength += alpha[row] * alpha[row];
        }
        alphaLength = sqrt(alphaLength);

        // True: We don't want to do the below since there isn't that amount of cols in Q and upper doesn't have col rows
        if (col >= rowAmount) continue; 
        
        upper[col][col] = alphaLength;

        // Calculatede the orthogonal column q
        for (int row = 0; row < alpha.size(); ++row) {
            orthogonal[row][col] = alpha[row] / alphaLength;
        }
    }
    return Matrix<T>::Qr(orthogonal, upper); 
}

// !!! unstabel might never converge
template <typename T>
typename Matrix<T>::Eigen Matrix<T>::calcEigen(double tolerance) const {
    if (!isSquare()) throw std::invalid_argument("Can not calculate Eigen value and vector of non square matrix");

    Matrix<T> A = Matrix<T>(*this);
    size_t rowAmount = A.getRowAmount();
    Matrix<T>::Qr qr = A.decompQR();
    Matrix<T> eigenVectors = qr.orthogonal; 

    std::vector<T> eigenValues;

    A = qr.upper * qr.orthogonal;

    while (!A.isUpperTriangular(tolerance)) {
        Matrix<T>::Qr qr = A.decompQR();
        A = qr.upper * qr.orthogonal ;
        eigenVectors = eigenVectors * qr.orthogonal;
    } 

    for (size_t i = 0; i < rowAmount; ++i) {
        eigenValues.push_back(A[i][i]);
    }
    return Matrix<T>::Eigen(eigenVectors, eigenValues);
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

                T currentElement = std::abs(upper[k][j]);
                if ((largestPivotRow == -1 && currentElement != 0) || 
                    (largestPivotRow != -1 && currentElement > std::abs(upper[largestPivotRow][j]))) {
                        largestPivotRow = k;
                }
            }
            pivotCol = j;
            if (largestPivotRow != -1) break; // no pivot only zeroes if this is the case
        }
        if (pivotCol > rowAmount) break;
        if (largestPivotRow == -1) continue;
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
    if (!isSquare()) throw std::runtime_error("Can not get inverse of non square matrix");

    size_t size = getRowAmount();
    Matrix<T> iden = getIdentity(size);
    Matrix<T> org(*this);

    for (int i = 0; i < size; ++i) {

        T scalar = org[i][i];

        bool isSingular = true;

        if (scalar == 0) {
            for (int j = i + 1; j < size; ++j) {
                if (org[j][i] != 0) {
                    isSingular = false;
                    org.swapRow(i, j);
                    iden.swapRow(i, j);
                    scalar = org[i][i];
                    break;
                }
            }
        }
        else isSingular = false;

        if (isSingular) throw  std::runtime_error("Cannot get inverse of singular matrix");
        for (int j = 0; j < size; ++j) { 
            org[i][j] = org[i][j] / scalar;
            iden[i][j] = iden[i][j] / scalar;
        }

        for (int j = 0; j < size; ++j) {
            if (i == j) continue;
            T scalar = org[j][i]; 
            if (scalar == 0) continue;
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
    if (getColAmount() != B.getRowAmount()) throw std::invalid_argument("Matrix dimensions do not allow multiplication");

    size_t rowAmount = getRowAmount();
    size_t colAmount = B.getColAmount();

    Matrix<T> product(rowAmount, colAmount, T()); 
    for (int row = 0; row < rowAmount; ++row) {
        for (int col = 0; col < colAmount; ++col) {
            T sum = T(); 
            for (int k = 0; k < colAmount; ++k) { 
                sum += data[row][k] * B[k][col];
            }
            product[row][col] = sum;
        }
    }
    return product;
}

template <typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& vec) const {
    if (getColAmount() != vec.size()) throw  std::invalid_argument("Matrix dimensions do not allow multiplication");
    
    size_t rowAmount = getRowAmount();

    std::vector<T> product(vec.size(), T()); 

    for (int i = 0; i < vec.size(); ++i) {
        T sum = T(); 
        for (int j = 0; j < rowAmount; ++j) {
            sum += data[i][j] * vec[j];
        }
        product[i] = sum;
    }

    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const  {
    size_t colAmount = getColAmount();
    size_t rowAmount = getRowAmount();

    Matrix<T> product(rowAmount, colAmount);
    for (size_t row = 0; row < rowAmount; ++row) {
        for (size_t col = 0; col < colAmount; ++col) {
            product[row][col] = scalar * data[row][col]; 
        }
    }
    return product;
}

// !!! doesn't work
template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, const Matrix<T>& B)  {
    Matrix<T> A;
    A.addRow(vec);
    Matrix<T> productMatrix = (A * B);
    std::vector<T> product = productMatrix.getCol(0);
    return product;
}

// !!! doesn't work
template <typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& B)  {
    return B * scalar;
}

template <typename T>
Matrix<T> Matrix<T>::pow(const T& exponent) const { // !!! Problem with getting correct eigen vectors 
    if (!isSquare()) throw std::invalid_argument("Can not take power of non square matrix");

    Matrix<T> A = Matrix<T>(*this);
    if (exponent < 0) A = A.inverse();
    T absoluteExponent = std::abs(exponent);
 
    Matrix<T>::Eigen eigen = A.calcEigen();
    std::vector<T> eigenValues = eigen.valueVec;
    size_t eigenValuesAmount = eigenValues.size();

    Matrix<T> eigenMatrix = eigen.vectorVec;
    Matrix<T> inverseEigenMatrix = eigenMatrix.inverse();
    Matrix<T> diagonal = Matrix<T>(eigenValuesAmount, eigenValuesAmount, static_cast<T>(0));

    for (int i = 0; i < eigenValuesAmount; ++i) {
        std::cout << eigenValues[i] <<  " "; 
    }

    for (int i = 0; i < eigenValuesAmount; ++i) {
        diagonal[i][i] = std::pow(eigenValues[i], absoluteExponent);
    }

    return (eigenMatrix * diagonal * inverseEigenMatrix);
}

// For non-const objects
template <typename T>
std::vector<T>& Matrix<T>::operator[](int index) {
    if (index >= rowAmount || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

// For const objects
template <typename T>
const std::vector<T>& Matrix<T>::operator[](int index) const {
    if (index >= rowAmount || index < 0) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

template <typename T>
std::vector<T> Matrix<T>::getCol(size_t colIndex) const {
    if (colIndex > colAmount -1) throw std::out_of_range("Invalid col index");
    std::vector<T> product;
    for (int i = 0; i < rowAmount; ++i) {
        product.push_back(data[i][colIndex]);
    } 
    return product;
}

template <typename T>
std::vector<T> Matrix<T>::getRow(size_t rowIndex) const {
    if (rowIndex > rowAmount -1) throw std::out_of_range("Invalid row index");
    return data[rowIndex];
}

template <typename T>
const T& Matrix<T>::getElement(size_t row, size_t col) const {
    if (row > rowAmount - 1 || col > colAmount - 1) throw std::out_of_range("Invalid row or col index");
    return data[row][col];
}

template <typename T>
T& Matrix<T>::getElement(size_t row, size_t col) {
    if (row > rowAmount - 1 || col > colAmount - 1) throw std::out_of_range("Invalid row or col index");
    return data[row][col];
}

template <typename T>
void Matrix<T>::assignElement(size_t row, size_t col, T el) {
    if (row > rowAmount - 1 || col > colAmount - 1) throw std::out_of_range("Invalid row or col index");
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
    if (rowA > rowAmount - 1 || rowB > rowAmount - 1) throw std::runtime_error("Tried to swap non existing row");
    std::vector<T> tempRow = data[rowA];
    (*this).setRow(rowA, data[rowB]);
    (*this).setRow(rowB, tempRow);
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::swapCol(size_t colA, size_t colB) {
    if (colA > colAmount - 1 || colB > colAmount - 1) throw std::runtime_error("Tried to swap non existing col");
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