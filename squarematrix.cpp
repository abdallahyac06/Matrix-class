#include "squarematrix.h"
#include <stdexcept>

SquareMatrix::SquareMatrix(int rows): Matrix(rows, rows) {}

SquareMatrix::SquareMatrix(const SquareMatrix &other): Matrix(other) {}

SquareMatrix::SquareMatrix(SquareMatrix &&other): Matrix(std::move(other)) {}

SquareMatrix::SquareMatrix(const Matrix &other): Matrix(other) {
    if (getRows() != getCols()) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(Matrix &&other): Matrix(std::move(other)) {
    if (getRows() != getCols()) {
        throw std::invalid_argument("Matrix is not square.");
    }
}

SquareMatrix::SquareMatrix(vector<double*> values, int rows): Matrix(values, rows, rows) {}

double SquareMatrix::trace() const {
    double result = 0;
    for (int i = 0; i < getRows(); ++i) {
        result += data[i][i];
    }
    return result;
}

double SquareMatrix::determinant() const {
    SquareMatrix copy(*this);
    int r = 0, c = 0;
    double pivot, result = 1;

    while (c < getCols()) {
        while (r < getRows() && !copy[r][c]) {
            ++r;
        }
        if (r == getRows()) {
            return 0.0;
        }
        if (r != c) {
            std::swap(copy[r], copy[c]);
            result = -result;
        }

        pivot = copy[c][c];
        for (int i = c + 1; i < getRows(); ++i) {
            if (copy[i][c]) {
                copy.addMultipleRow(i, c, -copy[i][c] / pivot);
            }
        }
        r = ++c;
    }

    for (int i = 0; i < getRows(); ++i) {
        result *= copy[i][i];
    }
    return result;
}

double SquareMatrix::determinantRecursive() const {
    if (getRows() == 1) {
        return data[0][0];
    }

    double result = 0;
    for (int i = 0; i < getRows(); ++i) {
        result += (i & 1 ? -1 : 1) * data[0][i] * this->operator()(0, i).determinantRecursive();
    }

    return result;
}

bool SquareMatrix::isLTriangular() const {
    for (int i = 0; i < getRows() - 1; ++i) {
        for (int j = i + 1; j < getCols(); ++j) {
            if (data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool SquareMatrix::isUTriangular() const {
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
    return isLTriangular() && isUTriangular();
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

const SquareMatrix SquareMatrix::operator()(int row, int col) const {
    if (row < 0 || row >= getRows()) {
        throw std::out_of_range("Row index out of bounds.");
    }

    if (col < 0 || col >= getCols()) {
        throw std::out_of_range("Column index out of bounds.");
    }
    
    SquareMatrix result(getRows() - 1);
    for (int i = 0; i < getRows() - 1; ++i) {
        for (int j = 0; j < getCols() - 1; ++j) {
            result[i][j] = data[i + (i >= row)][j + (j >= col)];
        }
    }
    
    return result;
}

const SquareMatrix SquareMatrix::adjoint() const {
    SquareMatrix result(getRows());
    for (int i = 0; i < getRows(); ++i) {
        for (int j = 0; j < getCols(); ++j) {
            result[j][i] = ((i ^ j) & 1 ? -1 : 1) * this->operator()(i, j).determinant();
        }
    }

    return result;
}

const SquareMatrix SquareMatrix::inverse() const {
    double det = determinant();
    if (!det) {
        throw std::logic_error("This matrix is not invertible.");
    }

    return SquareMatrix(adjoint() / det);
}

SquareMatrix &SquareMatrix::operator=(const SquareMatrix &other) {
    this->Matrix::operator=(other);
    return *this;
}

const SquareMatrix operator*(double scalar, const SquareMatrix &matrix) {
    return SquareMatrix(matrix * scalar);
}

std::ostream& operator<<(std::ostream &os, const SquareMatrix &matrix) {
    return os << Matrix(matrix);
}

std::istream& operator>>(std::istream &is, SquareMatrix &matrix) {
    for (double *row: matrix.data) {
        for (int i = 0; i < matrix.getCols(); ++i) {
            is >> row[i];
        }
    }
    
    return is;
}