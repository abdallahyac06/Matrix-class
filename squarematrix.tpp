#include "squarematrix.hpp"
#include <vector>
#include <iostream>

using std::ostream;
using std::istream;

using std::move;
using std::swap;

template <typename T>
SquareMatrix<T>::SquareMatrix(unsigned long size, const T& zero) : Matrix<T>(size, size, zero) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T>& other) : Matrix<T>(other) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(SquareMatrix<T>&& other) : Matrix<T>(move(other)) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(const Matrix<T>& other) : Matrix<T>(other) {
    if (this->ROWS != this->COLS) {
        throw MatrixException("Matrix is not square.");
    }
}

template <typename T>
SquareMatrix<T>::SquareMatrix(Matrix<T>&& other) : Matrix<T>(move(other)) {
    if (this->ROWS != this->COLS) {
        throw MatrixException("Matrix is not square.");
    }
}

template <typename T>
SquareMatrix<T>::SquareMatrix(const T* const* values, unsigned long size, const T& zero) : Matrix<T>(values, size, size, zero) {}

template <typename T>
T SquareMatrix<T>::determinantRecursive() {
    if (this->ROWS == 1) {
        return this->data[0][0];
    }

    unsigned long r = 0;
    while (r < this->ROWS && this->data[r][0] == this->ZERO) {
        ++r;
    }
    if (r == this->ROWS) {
        return this->ZERO;
    }

    for (unsigned long i = 0; i < this->ROWS; ++i) {
        if (i != r && this->data[i][0] != this->ZERO) {
            this->addMultipleRow(i, r, -this->data[i][0] / this->data[r][0]);
        }
    }

    return (r & 1 ? -this->data[r][0] : this->data[r][0]) * operator()(r, 0).determinantRecursive();
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::transpose() const {
    return SquareMatrix<T>(Matrix<T>::transpose());
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::ref() const {
    SquareMatrix<T> result(*this);
    unsigned long r = 0, c = 0;

    while (c < this->COLS) {
        while (r < this->ROWS && result[r][c] == this->ZERO) {
            ++r;
        }
        if (r == this->ROWS) {
            r = c++;
            continue;
        }

        swap(result[r], result[c]);
        for (unsigned long i = c + 1; i < this->ROWS; ++i) {
            if (result[i][c] != this->ZERO) {
                result.addMultipleRow(i, c, -result[i][c] / result[c][c], c + 1);
                result[i][c] = this->ZERO;
            }
        }
        r = ++c;
    }
    
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::rref() const {
    SquareMatrix<T> result(*this);
    unsigned long r = 0, c = 0;

    while (c < this->COLS) {
        while (r < this->ROWS && result[r][c] == this->ZERO) {
            ++r;
        }
        if (r == this->ROWS) {
            r = c++;
            continue;
        }

        swap(result[r], result[c]);
        result.divideRow(c, result[c][c], c);
        for (unsigned long i = 0; i < this->ROWS; ++i) {
            if (result[i][c] != this->ZERO && i != c) {
                result.addMultipleRow(i, c, -result[i][c]);
            }
        }
        r = ++c;
    }

    return result;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator=(const SquareMatrix<T>& other) {
    Matrix<T>::operator=(other);
    return *this;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator+(const Matrix<T>& other) const {
    return SquareMatrix<T>(Matrix<T>::operator+(other));
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator-(const Matrix<T>& other) const {
    return SquareMatrix<T>(Matrix<T>::operator-(other));
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator-() const {
    return SquareMatrix<T>(Matrix<T>::operator-());
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const SquareMatrix<T>& other) const {
    return SquareMatrix<T>(Matrix<T>::operator*(other));
}

template <typename T>
Matrix<T> SquareMatrix<T>::operator*(const Matrix<T>& other) const {
    return Matrix<T>::operator*(other);
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(T scalar) const {
    return SquareMatrix<T>(Matrix<T>::operator*(scalar));
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator/(T scalar) const {
    return SquareMatrix<T>(Matrix<T>::operator/(scalar));
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator+=(const Matrix<T>& other) {
    return operator=(operator+(other));
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator-=(const Matrix<T>& other) {
    return operator=(operator-(other));
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const Matrix<T>& other) {
    return operator=(operator*(other));
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(T scalar) {
    return operator=(operator*(scalar));
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator/=(T scalar) {
    return operator=(operator/(scalar));
}

template <typename T>
const SquareMatrix<T> SquareMatrix<T>::id(unsigned long size, T scalar) {
    SquareMatrix<T> result(size);
    for (unsigned long i = 0; i < size; ++i) {
        result[i][i] = scalar;
    }

    return result;
}

template <typename T>
T SquareMatrix<T>::trace() const {
    T result = this->ZERO;
    for (unsigned long i = 0; i < this->ROWS; ++i) {
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
    for (unsigned long i = 0; i < this->ROWS - 1; ++i) {
        for (unsigned long j = i + 1; j < this->COLS; ++j) {
            if (this->data[i][j] != this->ZERO) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
bool SquareMatrix<T>::isUTriangular() const {
    for (unsigned long i = 1; i < this->ROWS; ++i) {
        for (unsigned long j = 0; j < i; ++j) {
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
    for (unsigned long i = 1; i < this->ROWS; ++i) {
        for (unsigned long j = 0; j < i; ++j) {
            if (this->data[i][j] != this->data[j][i]) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator()(unsigned long row, unsigned long col) const {
    if (row < 0 || row >= this->ROWS) {
        throw MatrixException("Row index out of bounds.");
    }

    if (col < 0 || col >= this->COLS) {
        throw MatrixException("Column index out of bounds.");
    }
    
    SquareMatrix<T> result(this->ROWS - 1);
    for (unsigned long i = 0; i < this->ROWS - 1; ++i) {
        for (unsigned long j = 0; j < this->COLS - 1; ++j) {
            result[i][j] = this->data[i + (i >= row)][j + (j >= col)];
        }
    }
    
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::adjoint() const {
    SquareMatrix<T> result(this->ROWS);
    for (unsigned long i = 0; i < this->ROWS; ++i) {
        for (unsigned long j = 0; j < this->COLS; ++j) {
            result[j][i] = operator()(i, j).determinantRecursive();
            if ((i ^ j) & 1) {
                result[j][i] = -result[j][i];
            }
        }
    }

    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::inverse() const {
    SquareMatrix<T> copy(*this), result(id(this->ROWS));
    unsigned long r = 0, c = 0;

    while (c < this->COLS) {
        while (r < this->ROWS && copy[r][c] == this->ZERO) {
            ++r;
        }
        if (r == this->ROWS) {
            throw MatrixException("Matrix is singular.");
        }

        swap(result.data[r], result.data[c]);
        swap(copy.data[r], copy.data[c]);
        result.divideRow(c, copy[c][c]);
        copy.divideRow(c, copy[c][c]);

        for (unsigned long i = 0; i < this->ROWS; ++i) {
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
ostream& operator<<(ostream& os, const SquareMatrix<T>& matrix) {
    return os << Matrix<T>(matrix);
}

template <typename T>
istream& operator>>(istream& is, SquareMatrix<T>& matrix) {
    return is >> static_cast<Matrix<T>&>(matrix);
}
    
template <typename T>
SquareMatrix<T> operator*(T scalar, const SquareMatrix<T>& matrix) {
    return matrix * scalar;
}