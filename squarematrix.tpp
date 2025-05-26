#include "squarematrix.hpp"
#include <stdexcept>

template <typename T>
SquareMatrix<T>::SquareMatrix(int size, const T &zero): Matrix<T>(size, size, zero) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T> &other): Matrix<T>(other) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(SquareMatrix<T> &&other): Matrix<T>(std::move(other)) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(const Matrix<T> &other): Matrix<T>(other) {
    if (this->getRows() != this->getCols()) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

template <typename T>
SquareMatrix<T>::SquareMatrix(Matrix<T> &&other): Matrix<T>(std::move(other)) {
    if (this->getRows() != this->getCols()) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

template <typename T>
SquareMatrix<T>::SquareMatrix(vector<T*> values, int size, const T &zero): Matrix<T>(values, size, size, zero) {}

template <typename T>
T SquareMatrix<T>::determinantRecursive() {
    if (this->getRows() == 1) {
        return this->data[0][0];
    }

    int r = 0;
    while (r < this->getRows() && this->data[r][0] == this->ZERO) {
        ++r;
    }
    if (r == this->getRows()) {
        return this->ZERO;
    }

    for (int i = 0; i < this->getRows(); ++i) {
        if (i != r && this->data[i][0] != this->ZERO) {
            this->addMultipleRow(i, r, -this->data[i][0] / this->data[r][0]);
        }
    }

    return (r & 1 ? -this->data[r][0] : this->data[r][0]) * operator()(r, 0).determinantRecursive();
}

template <typename T>
const SquareMatrix<T> SquareMatrix<T>::id(int size) {
    SquareMatrix<T> result(size);
    for (int i = 0; i < size; ++i) {
        result[i][i] = T(1);
    }

    return result;
}

template <typename T>
T SquareMatrix<T>::trace() const {
    T result = this->ZERO;
    for (int i = 0; i < this->getRows(); ++i) {
        result += this->data[i][i];
    }

    return result;
}

template <typename T>
T SquareMatrix<T>::determinant() const {
    return SquareMatrix<T>(*this).determinantRecursive();
}

template <typename T>
bool SquareMatrix<T>::isLTriangular() const {
    for (int i = 0; i < this->getRows() - 1; ++i) {
        for (int j = i + 1; j < this->getCols(); ++j) {
            if (this->data[i][j] != this->ZERO) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
bool SquareMatrix<T>::isUTriangular() const {
    for (int i = 1; i < this->getRows(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (this->data[i][j] != this->ZERO) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
bool SquareMatrix<T>::isDiagonal() const {
    return isLTriangular() && isUTriangular();
}

template <typename T>
bool SquareMatrix<T>::isSymmetric() const {
    for (int i = 1; i < this->getRows(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (this->data[i][j] != this->data[j][i]) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator()(int row, int col) const {
    if (row < 0 || row >= this->getRows()) {
        throw std::out_of_range("Row index out of bounds.");
    }

    if (col < 0 || col >= this->getCols()) {
        throw std::out_of_range("Column index out of bounds.");
    }
    
    SquareMatrix<T> result(this->getRows() - 1);
    for (int i = 0; i < this->getRows() - 1; ++i) {
        for (int j = 0; j < this->getCols() - 1; ++j) {
            result[i][j] = this->data[i + (i >= row)][j + (j >= col)];
        }
    }
    
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::adjoint() const {
    SquareMatrix<T> result(this->getRows());
    for (int i = 0; i < this->getRows(); ++i) {
        for (int j = 0; j < this->getCols(); ++j) {
            result[j][i] = operator()(i, j).determinant();
            if ((i ^ j) & 1) {
                result[j][i] = -result[j][i];
            }
        }
    }

    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::inverse() const {
    SquareMatrix<T> copy(*this), result(id(this->getRows()));
    int r = 0, c = 0;

    while (c < this->getCols()) {
        while (r < this->getRows() && copy[r][c] == this->ZERO) {
            ++r;
        }
        if (r == this->getRows()) {
            throw std::logic_error("Matrix is singular.");
        }

        std::swap(result.data[r], result.data[c]);
        std::swap(copy.data[r], copy.data[c]);
        result.divideRow(c, copy[c][c]);
        copy.divideRow(c, copy[c][c]);

        for (int i = 0; i < this->getRows(); ++i) {
            if (copy[i][c] != this->ZERO && i != c) {
                result.addMultipleRow(i, c, -copy[i][c]);
                copy.addMultipleRow(i, c, -copy[i][c]);
            }
        }
        r = ++c;
    }

    return result;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator=(const SquareMatrix<T> &other) {
    Matrix<T>::operator=(other);
    return *this;
}

template <typename T>
const SquareMatrix<T> operator*(T scalar, const SquareMatrix<T> &matrix) {
    return SquareMatrix<T>(matrix * scalar);
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const SquareMatrix<T> &matrix) {
    return os << Matrix<T>(matrix);
}

template <typename T>
std::istream& operator>>(std::istream &is, SquareMatrix<T> &matrix) {
    for (T *row: matrix.data) {
        for (int i = 0; i < matrix.getCols(); ++i) {
            is >> row[i];
        }
    }
    
    return is;
}