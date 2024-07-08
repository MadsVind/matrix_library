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
    if (!isSameSize(B)) {
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
    if (!isSameSize(B)) {
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
Matrix<T> Matrix<T>::ref(double tolerance) const {
    Matrix<T> ref(*this);

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    
    T zero = static_cast<T>(0);
    int pivotCol = 0;
    for (int i = 0; i < rowAmount; ++i) { // Find pivot in every row
        int largestPivotRow = -1;
        for (int j = pivotCol; j < colAmount; ++j) { // Go to next column if column has no pivot
            for (int k = i; k < rowAmount; ++k) { 
                T currentElement = std::abs(ref[k][j]);
                if ((largestPivotRow == -1 && currentElement > tolerance)) {    
                    largestPivotRow = k;
                }
                if (largestPivotRow != -1 && currentElement > std::abs(ref[largestPivotRow][j])) largestPivotRow = k;
            }
            pivotCol = j;
            if (largestPivotRow != -1) break; // no pivot only zeroes if this is the case
        }

        if (pivotCol > rowAmount) break;
        if (largestPivotRow == -1) continue;
        if (largestPivotRow != i) ref.swapRow(i, largestPivotRow);
        T scalar = ref[i][pivotCol];

        for (int j = 0; j < colAmount; ++j) { 
            ref[i][j] = ref[i][j] / scalar;
        }
        

        for (int j = i + 1; j < rowAmount; ++j) {
            T scalar = ref[j][i]; 
            for (int k = 0; k < colAmount; ++k) {
                ref[j][k] -= scalar * ref[i][k];
            }
        }
    }
    return ref;
}

template <typename T>
Matrix<T> Matrix<T>::rref(double tolerance) const {
    Matrix<T> rref((*this).ref());

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();
    
    int pivotCol = 0;
    for (int i = 0; i < rowAmount; ++i) { // Find pivot in every row
        for (int j = 0; j < rowAmount; ++j) {
            if (i == j) continue;
            T scalar = rref[j][i]; 
            for (int k = i; k < colAmount; ++k) {
                rref[j][k] -= scalar * rref[i][k];
            }
        }
    }
    return rref;
}

template <typename T>
std::vector<T> backSubstitution(Matrix<T> A, double tolerance) { // !!! does not work fix
    size_t rowAmount = A.getRowAmount();
    std::vector<T> res(rowAmount);
    for (int i = rowAmount - 1; 0 <= i; --i) {
        if (std::abs(A[i][i]) < tolerance) {
            res[i] = static_cast<T>(1); 
            continue;
        }
        int rowAmountIndex = (rowAmount - 1);
        int varAmount =  rowAmountIndex - i;
        T varSum = static_cast<T>(0);
        for (int j = 0; j < varAmount; ++j) {
            varSum += A[i][ rowAmountIndex - j] * res[rowAmountIndex - j];
        }
        res[i] = -varSum;
    }
    return res;
}

template <typename T>
typename Matrix<T>::Eigen Matrix<T>::calcEigen(double tolerance) const {
    tolerance * 10;
    if (!isSquare()) {
        std::cerr << "Can not calculate Eigen value and vector of non square matrix\n";
        return Matrix<T>::Eigen(std::vector<std::vector<T>>(), std::vector<T>());
    }
    Matrix<T> A = Matrix<T>(*this);
    do {
        Matrix<T>::Qr qr = A.decompQR();

        A = qr.ortogonal * qr.upper;

    } while (!A.isUpperTriangular(tolerance * 0.1)); // nessesary for tolerance to work kinda

    std::vector<T> eigenValues;
    size_t colAmount = A.getColAmount();
    size_t rowAmount = A.getRowAmount();

    for (size_t i = 0; i < rowAmount; ++i) {
        T rounded = std::round(A[i][i] / tolerance);
        eigenValues.push_back(rounded * tolerance);
    }
    A = Matrix<T>(*this); // think this is the thing which needed to be done back to drawing board
 
    Matrix<T> I = getIdentity(rowAmount);
    std::vector<std::vector<T>> eigenVectors;
    size_t valueAmount = eigenValues.size(); // probably should just be row amount
    for (size_t i = 0; i < valueAmount; ++i) {
        Matrix<T> linMatrix = A - (I * eigenValues[i]);
        linMatrix = linMatrix.ref(tolerance);
        eigenVectors.push_back(backSubstitution(linMatrix, tolerance));
    }

    return Matrix<T>::Eigen(eigenVectors, eigenValues); // !!!remember to change
}

template <typename T>
typename Matrix<T>::Qr Matrix<T>::decompQR() const {
    Matrix<T> orthogonal = Matrix<T>(rowAmount, rowAmount);
    Matrix<T> upper = Matrix<T>(rowAmount, colAmount, static_cast<T>(0));

    size_t rowAmount = getRowAmount();
    size_t colAmount = getColAmount();

    for (int i = 0; i < rowAmount; ++i) {
        std::vector<T> a = (*this).getCol(i);
        std::vector<T> scalars;

        // Calculate the scalars of the new orthogonal vector a (alpha)
        // These scalars are also the non digonal entries in upper
        for (int j = 0; j < i; ++j) {
            T scalar = static_cast<T>(0);
            std::vector<T> q = orthogonal.getCol(j);

            for (int k = 0; k < a.size(); ++k) {
                scalar += a[k] * q[k];
            }
            scalars.push_back(scalar);
            upper[j][i] = scalar;
            
        }

        // calculatede the new orthogonal vector a (alpha)
        std::vector<T> alpha;
        for (int j = 0; j < a.size(); ++j) {
            T el = a[j];
            for (int k = 0; k < i; ++k) {
                std::vector<T> q = orthogonal.getCol(k);
                el -= scalars[k] * q[j];
            }
            alpha.push_back(el);
        }
        
        //  calculate the length of new orthogonal vector a (alpha)
        T alphaLength = static_cast<T>(0);
        for (int j = 0; j < alpha.size(); ++j) {
            alphaLength += alpha[j] * alpha[j];
        }
        alphaLength = sqrt(alphaLength);
        upper[i][i] = alphaLength;

        // calculatede the orthogonal column q
        for (int j = 0; j < alpha.size(); ++j) {
            orthogonal[j][i] = alpha[j] / alphaLength;
        }
    }
    return Matrix<T>::Qr(orthogonal, upper);
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
    if (!isSquare()) {
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