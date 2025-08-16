#include "squarematrix.h"

SquareMatrix::SquareMatrix(int size): Matrix(size, size) {}

SquareMatrix::SquareMatrix(const SquareMatrix &other): Matrix(other) {}

SquareMatrix::SquareMatrix(SquareMatrix &&other): Matrix(std::move(other)) {}

SquareMatrix::SquareMatrix(const Matrix &other): Matrix(other) {
    if (getRows() != getCols()) {
        throw MatrixException("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(Matrix &&other): Matrix(std::move(other)) {
    if (getRows() != getCols()) {
        throw MatrixException("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(vector<double *> values, int size): Matrix(values, size, size) {}

double SquareMatrix::determinantRecursive() {
    if (getRows() == 1) {
        return data[0][0];
    }

    int r = 0;
    while (r < getRows() && !data[r][0]) {
        ++r;
    }
    if (r == getRows()) {
        return 0.0;
    }

    for (int i = 0; i < getRows(); ++i) {
        if (i != r && data[i][0]) {
            addMultipleRow(i, r, -data[i][0] / data[r][0]);
        }
    }

    return (r & 1 ? -data[r][0]: data[r][0]) * operator()(r, 0).determinantRecursive();
}

SquareMatrix SquareMatrix::id(int size) {
    SquareMatrix result(size);
    for (int i = 0; i < size; ++i) {
        result[i][i] = 1.0;
    }

    return result;
}

double SquareMatrix::trace() const {
    double result = 0;
    for (int i = 0; i < getRows(); ++i) {
        result += data[i][i];
    }

    return result;
}

double SquareMatrix::determinant() const {
    return SquareMatrix(*this).determinantRecursive();
}

bool SquareMatrix::isLowerTriangular() const {
    for (int i = 0; i < getRows() - 1; ++i) {
        for (int j = i + 1; j < getCols(); ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isUpperTriangular() const {
    for (int i = 1; i < getRows(); ++i) {
        for (int j = 0; j < i; ++j) {
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
    for (int i = 1; i < getRows(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (data[i][j] != data[j][i]) {
                return false;
            }
        }
    }

    return true;
}

SquareMatrix SquareMatrix::operator()(int row, int col) const {
    if (row < 0 || row >= getRows()) {
        throw MatrixException("Row index out of bounds.");
    }

    if (col < 0 || col >= getCols()) {
        throw MatrixException("Column index out of bounds.");
    }

    SquareMatrix result(getRows() - 1);
    for (int i = 0; i < getRows() - 1; ++i) {
        for (int j = 0; j < getCols() - 1; ++j) {
            result[i][j] = data[i + (i >= row)][j + (j >= col)];
        }
    }

    return result;
}

SquareMatrix SquareMatrix::adjoint() const {
    SquareMatrix result(getRows());
    for (int i = 0; i < getRows(); ++i) {
        for (int j = 0; j < getCols(); ++j) {
            result[j][i] = ((i ^ j) & 1 ? -1.0 : 1.0) * operator()(i, j).determinant();
        }
    }

    return result;
}

SquareMatrix SquareMatrix::inverse() const {
    SquareMatrix copy(*this), result(id(getRows()));
    int r = 0, c = 0;

    while (c < getCols()) {
        while (r < getRows() && !copy[r][c]) {
            ++r;
        }
        if (r == getRows()) {
            throw MatrixException("Matrix is singular.");
        }

        std::swap(result.data[r], result.data[c]);
        std::swap(copy.data[r], copy.data[c]);
        result.divideRow(c, copy[c][c]);
        copy.divideRow(c, copy[c][c]);

        for (int i = 0; i < getRows(); ++i) {
            if (copy[i][c] && i != c) {
                result.addMultipleRow(i, c, -copy[i][c]);
                copy.addMultipleRow(i, c, -copy[i][c]);
            }
        }
        r = ++c;
    }

    return result;
}

SquareMatrix &SquareMatrix::operator=(const SquareMatrix &other) {
    Matrix::operator=(other);
    return *this;
}

SquareMatrix operator*(double scalar, const SquareMatrix &matrix) {
    return SquareMatrix(matrix * scalar);
}

std::ostream &operator<<(std::ostream &os, const SquareMatrix &matrix) {
    return os << Matrix(matrix);
}

std::istream &operator>>(std::istream &is, SquareMatrix &matrix) {
    for (double *row: matrix.data) {
        for (int i = 0; i < matrix.getCols(); ++i) {
            is >> row[i];
        }
    }

    return is;
}