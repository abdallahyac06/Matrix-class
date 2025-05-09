#include "matrix.h"
#include <iostream>
#include <vector>
#include <stdexcept>

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

template <typename T>
class SquareMatrix : public Matrix<T> {
    public:
        SquareMatrix(int size = 4, const T &zero = T());
        SquareMatrix(const SquareMatrix<T> &other);
        SquareMatrix(SquareMatrix<T> &&other);
        SquareMatrix(const Matrix<T> &other);
        SquareMatrix(Matrix<T> &&other);
        SquareMatrix(vector<T*> values, int size, const T &zero = T());
        virtual ~SquareMatrix() = default;

        T trace() const;
        T determinant() const;
        T determinantRecursive() const;
        bool isLTriangular() const;
        bool isUTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        const SquareMatrix<T> operator()(int row, int col) const;
        const SquareMatrix<T> adjoint() const;
        const SquareMatrix<T> inverse() const;
        SquareMatrix<T> &operator=(const SquareMatrix<T> &other);

    template <typename U>
    friend const SquareMatrix<U> operator*(U scalar, const SquareMatrix<U> &matrix);
    template <typename U>
    friend ostream &operator<<(ostream &os, const SquareMatrix<U> &matrix);
    template <typename U>
    friend istream &operator>>(istream &is, SquareMatrix<U> &matrix);
};


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
T SquareMatrix<T>::trace() const {
    T result = this->ZERO;
    for (int i = 0; i < this->getRows(); ++i) {
        result += this->data[i][i];
    }
    return result;
}

template <typename T>
T SquareMatrix<T>::determinant() const {
    SquareMatrix<T> copy(*this);
    int r = 0, c = 0;
    T pivot;
    bool even_swaps = true;

    while (c < this->getCols()) {
        while (r < this->getRows() && copy[r][c] == this->ZERO) {
            ++r;
        }
        if (r == this->getRows()) {
            return this->ZERO;
        }
        if (r != c) {
            std::swap(copy[r], copy[c]);
            even_swaps = !even_swaps;
        }

        pivot = copy[c][c];
        for (int i = c + 1; i < this->getRows(); ++i) {
            if (copy[i][c] != this->ZERO) {
                copy.addMultipleRow(i, c, -copy[i][c] / pivot);
            }
        }
        r = ++c;
    }

    T result = even_swaps ? copy[0][0] : -copy[0][0];
    for (int i = 1; i < this->getRows(); ++i) {
        result *= copy[i][i];
    }

    return result;
}

template <typename T>
T SquareMatrix<T>::determinantRecursive() const {
    if (this->getRows() == 1) {
        return this->data[0][0];
    }

    T result = this->ZERO;
    for (int i = 0; i < this->getRows(); i += 2) {
        result += this->data[0][i] * operator()(0, i).determinantRecursive();
    }
    for (int i = 1; i < this->getRows(); i += 2) {
        result -= this->data[0][i] * operator()(0, i).determinantRecursive();
    }

    return result;
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
const SquareMatrix<T> SquareMatrix<T>::operator()(int row, int col) const {
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
const SquareMatrix<T> SquareMatrix<T>::adjoint() const {
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
const SquareMatrix<T> SquareMatrix<T>::inverse() const {
    T det = determinant();
    if (det == this->ZERO) {
        throw std::logic_error("This matrix is not invertible.");
    }

    return SquareMatrix<T>(adjoint() / det);
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator=(const SquareMatrix<T> &other) {
    this->Matrix<T>::operator=(other);
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

#endif