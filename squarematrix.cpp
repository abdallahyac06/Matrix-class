#include "squarematrix.h"

using std::move;
using std::swap;


SquareMatrix::SquareMatrix(unsigned long size): Matrix(size, size) {}

SquareMatrix::SquareMatrix(const SquareMatrix& other): Matrix(other) {}

SquareMatrix::SquareMatrix(SquareMatrix&& other): Matrix(move(other)) {}

SquareMatrix::SquareMatrix(const Matrix& other): Matrix(other) {
    if (getRows() != getCols()) {
        throw MatrixException("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(Matrix&& other): Matrix(move(other)) {
    if (getRows() != getCols()) {
        throw MatrixException("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(const double* const* values, unsigned long size): Matrix(values, size, size) {}

double SquareMatrix::determinantRecursive() {
    if (getRows() == 1) {
        return data[0][0];
    }

    unsigned long r = 0;
    while (r < getRows() && !data[r][0]) {
        ++r;
    }
    if (r == getRows()) {
        return 0.0;
    }

    for (unsigned long i = 0; i < getRows(); ++i) {
        if (i != r && data[i][0]) {
            addMultipleRow(i, r, -data[i][0] / data[r][0]);
        }
    }

    return (r & 1 ? -data[r][0]: data[r][0]) * operator()(r, 0).determinantRecursive();
}

SquareMatrix SquareMatrix::transpose() const {
    return SquareMatrix(Matrix::transpose());
}

SquareMatrix SquareMatrix::ref() const {
    SquareMatrix result(*this);
    unsigned long r = 0, c = 0;

    while (c < getCols()) {
        while (r < getRows() && !result[r][c]) {
            ++r;
        }
        if (r == getRows()) {
            r = c++;
            continue;
        }

        swap(result[r], result[c]);
        for (unsigned long i = c + 1; i < getRows(); ++i) {
            if (result[i][c]) {
                result.addMultipleRow(i, c, -result[i][c] / result[c][c], c + 1);
                result[i][c] = 0;
            }
        }
        r = ++c;
    }
    
    return result;
}

SquareMatrix SquareMatrix::rref() const {
    SquareMatrix result(*this);
    unsigned long r = 0, c = 0;

    while (c < getCols()) {
        while (r < getRows() && !result[r][c]) {
            ++r;
        }
        if (r == getRows()) {
            r = c++;
            continue;
        }

        swap(result[r], result[c]);
        result.divideRow(c, result[c][c], c);
        for (unsigned long i = 0; i < getRows(); ++i) {
            if (result[i][c] && i != c) {
                result.addMultipleRow(i, c, -result[i][c]);
            }
        }
        r = ++c;
    }

    return result;
}

SquareMatrix& SquareMatrix::operator=(const SquareMatrix& other) {
    Matrix::operator=(other);
    return *this;
}

SquareMatrix SquareMatrix::operator+(const Matrix& other) const {
    return SquareMatrix(Matrix::operator+(other));
}

SquareMatrix SquareMatrix::operator-(const Matrix& other) const {
    return SquareMatrix(Matrix::operator-(other));
}

SquareMatrix SquareMatrix::operator-() const {
    return SquareMatrix(Matrix::operator-());
}

SquareMatrix SquareMatrix::operator*(const SquareMatrix& other) const {
    return SquareMatrix(Matrix::operator*(other));
}

Matrix SquareMatrix::operator*(const Matrix& other) const {
    return Matrix::operator*(other);
}

SquareMatrix SquareMatrix::operator*(double scalar) const {
    return SquareMatrix(Matrix::operator*(scalar));
}

SquareMatrix SquareMatrix::operator/(double scalar) const {
    return SquareMatrix(Matrix::operator/(scalar));
}

SquareMatrix& SquareMatrix::operator+=(const Matrix& other) {
    return operator=(operator+(other));
}

SquareMatrix& SquareMatrix::operator-=(const Matrix& other) {
    return operator=(operator-(other));
}

SquareMatrix& SquareMatrix::operator*=(const Matrix& other) {
    return operator=(operator*(other));
}

SquareMatrix& SquareMatrix::operator*=(double scalar) {
    return operator=(operator*(scalar));
}

SquareMatrix& SquareMatrix::operator/=(double scalar) {
    return operator=(operator/(scalar));
}

SquareMatrix SquareMatrix::id(unsigned long size) {
    SquareMatrix result(size);
    for (unsigned long i = 0; i < size; ++i) {
        result[i][i] = 1.0;
    }

    return result;
}

double SquareMatrix::trace() const {
    double result = 0;
    for (unsigned long i = 0; i < getRows(); ++i) {
        result += data[i][i];
    }

    return result;
}

double SquareMatrix::determinant() const {
    return SquareMatrix(*this).determinantRecursive();
}

bool SquareMatrix::isLowerTriangular() const {
    for (unsigned long i = 0; i < getRows() - 1; ++i) {
        for (unsigned long j = i + 1; j < getCols(); ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isUpperTriangular() const {
    for (unsigned long i = 1; i < getRows(); ++i) {
        for (unsigned long j = 0; j < i; ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isDiagonal() const {
    return isLowerTriangular() && isUpperTriangular();
}

bool SquareMatrix::isSymmetric() const {
    for (unsigned long i = 1; i < getRows(); ++i) {
        for (unsigned long j = 0; j < i; ++j) {
            if (data[i][j] != data[j][i]) {
                return false;
            }
        }
    }

    return true;
}

SquareMatrix SquareMatrix::operator()(unsigned long row, unsigned long col) const {
    if (row < 0 || row >= getRows()) {
        throw MatrixException("Row index out of bounds.");
    }

    if (col < 0 || col >= getCols()) {
        throw MatrixException("Column index out of bounds.");
    }

    SquareMatrix result(getRows() - 1);
    for (unsigned long i = 0; i < getRows() - 1; ++i) {
        for (unsigned long j = 0; j < getCols() - 1; ++j) {
            result[i][j] = data[i + (i >= row)][j + (j >= col)];
        }
    }

    return result;
}

SquareMatrix SquareMatrix::adjoint() const {
    SquareMatrix result(getRows());
    for (unsigned long i = 0; i < getRows(); ++i) {
        for (unsigned long j = 0; j < getCols(); ++j) {
            result[j][i] = ((i ^ j) & 1 ? -1.0 : 1.0) * operator()(i, j).determinant();
        }
    }

    return result;
}

SquareMatrix SquareMatrix::inverse() const {
    SquareMatrix copy(*this), result(id(getRows()));
    unsigned long r = 0, c = 0;

    while (c < getCols()) {
        while (r < getRows() && !copy[r][c]) {
            ++r;
        }
        if (r == getRows()) {
            throw MatrixException("Matrix is singular.");
        }

        swap(result[r], result[c]);
        swap(copy[r], copy[c]);
        result.divideRow(c, copy[c][c]);
        copy.divideRow(c, copy[c][c]);

        for (unsigned long i = 0; i < getRows(); ++i) {
            if (copy[i][c] && i != c) {
                result.addMultipleRow(i, c, -copy[i][c]);
                copy.addMultipleRow(i, c, -copy[i][c]);
            }
        }
        r = ++c;
    }

    return result;
}

SquareMatrix operator*(double scalar, const SquareMatrix& matrix) {
    return SquareMatrix(matrix * scalar);
}

ostream& operator<<(ostream& os, const SquareMatrix& matrix) {
    return os << static_cast<Matrix>(matrix);
}

istream& operator>>(istream& is, SquareMatrix& matrix) {
    return is >> static_cast<Matrix&>(matrix);
}