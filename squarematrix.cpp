#include "squarematrix.h"
#include <iostream>
#include <stdexcept>

SquareMatrix::SquareMatrix(int rows): Matrix(rows, rows) {}

SquareMatrix::SquareMatrix(const SquareMatrix &other): Matrix(other.data, other.ROWS, other.ROWS) {}

SquareMatrix::SquareMatrix(SquareMatrix &&other): Matrix(std::move(other)) {}

SquareMatrix::SquareMatrix(const Matrix &other): SquareMatrix(other.getData(), other.getRows()) {
    if (ROWS != COLS) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(Matrix &&other): Matrix(std::move(other)) {
    if (ROWS != COLS) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(double **array, int rows): Matrix(array, rows, rows) {}

double SquareMatrix::trace() const {
    double result = 0;
    for (int i = 0; i < ROWS; ++i) {
        result += data[i][i];
    }
    return result;
}

double SquareMatrix::determinant() const {
    SquareMatrix copy(*this);
    int r = 0, c = 0;
    double pivot, result = 1;

    while (c < COLS) {
        while (r < ROWS && !copy[r][c]) {
            ++r;
        }
        if (r == ROWS) {
            return 0.0;
        }
        if (r != c) {
            copy.swapRows(r, c);
            result = -result;
        }

        pivot = copy[c][c];
        for (int i = c + 1; i < ROWS; ++i) {
            if (copy[i][c]) {
                copy.addMultipleRow(i, c, -copy[i][c] / pivot);
            }
        }
        r = ++c;
    }

    for (int i = 0; i < ROWS; ++i) {
        result *= copy[i][i];
    }
    return result;
}

double SquareMatrix::determinantRecursive() const {
    if (ROWS == 1) {
        return data[0][0];
    }

    double result = 0;
    for (int i = 0; i < ROWS; ++i) {
        result += ((i & 1) ? -1 : 1) * data[0][i] * SquareMatrix(this->operator()(0, i)).determinantRecursive();
    }

    return result;
}

bool SquareMatrix::isLTriangular() const {
    for (int i = 0; i < ROWS - 1; ++i) {
        for (int j = i + 1; j < COLS; ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isUTriangular() const {
    for (int i = 1; i < ROWS; ++i) {
        for (int j = 0; j < i; ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isDiagonal() const {
    return isLTriangular() && isUTriangular();
}

bool SquareMatrix::isSymmetric() const {
    for (int i = 1; i < ROWS; ++i) {
        for (int j = 0; j < i; ++j) {
            if (data[i][j] != data[j][i]) {
                return false;
            }
        }
    }

    return true;
}

SquareMatrix SquareMatrix::adjoint() const {
    SquareMatrix result(ROWS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[j][i] = ((i & j & 1) ? -1 : 1) * SquareMatrix(this->operator()(i, j)).determinant();
        }
    }

    return result;
}

SquareMatrix SquareMatrix::inverse() const {
    double det = determinant();
    if (!det) {
        throw std::logic_error("This matrix is not invertible.");
    }

    return 1 / det * adjoint();
}

SquareMatrix operator*(double scalar, const SquareMatrix &matrix) {
    SquareMatrix result(matrix);
    for (int i = 0; i < result.ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }
    
    return result;
}

std::ostream& operator<<(std::ostream &os, const SquareMatrix &matrix) {
    return os << Matrix(matrix);
}

std::istream& operator>>(std::istream &is, SquareMatrix &matrix) {
    for (int i = 0; i < matrix.ROWS; ++i) {
        for (int j = 0; j < matrix.COLS; ++j) {
            is >> matrix[i][j];
        }
    }
    
    return is;
}