#include "matrix.h"
#include <iostream>
#include <stdexcept>

Matrix::Matrix(int rows, int cols): ROWS(rows), COLS(cols) {
    if (ROWS < 1 || COLS < 1) {
        throw std::invalid_argument("Matrix dimensions must be greater than 0.");
    }

    data = new double*[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        data[i] = new double[COLS];
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix &other): Matrix(other.ROWS, other.COLS) {
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = other[i][j];
        }
    }
}

Matrix::Matrix(const double **array, int rows, int cols): Matrix(rows, cols) {
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = array[i][j];
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < ROWS; ++i) {
        delete[] data[i];
        data[i] = nullptr;
    }
    delete[] data;
    data = nullptr;
}

void Matrix::checkArgs(int row, int col) const {
    if (row < 0 || row >= ROWS) {
        throw std::invalid_argument("Row index out of bounds.");
    }
    if (col < 0 || col >= COLS) {
        throw std::invalid_argument("Col index out of bounds.");
    }
}

void Matrix::setRow(int row, const double *values) {
    checkArgs(row, 0);

    for (int i = 0; i < COLS; ++i) {
        data[row][i] = values[i];
    }
}

void Matrix::setCol(int col, const double *values) {
    checkArgs(0, col);
    
    for (int i = 0; i < ROWS; ++i) {
        data[i][col] = values[i];
    }
}

void Matrix::swapRows(int row1, int row2) {
    checkArgs(row1, 0);
    checkArgs(row2, 0);

    for (int i = 0; i < COLS; ++i) {
        std::swap(data[row1][i], data[row2][i]);
    }
}

void Matrix::swapCols(int col1, int col2) {
    checkArgs(0, col1);
    checkArgs(0, col2);

    for (int i = 0; i < ROWS; ++i) {
        std::swap(data[i][col1], data[i][col2]);
    }
}

void Matrix::multiplyRow(int row, double scalar) {
    checkArgs(row, 0);

    for (int i = 0; i < COLS; ++i) {
        data[row][i] *= scalar;
    }
}

void Matrix::multiplyCol(int col, double scalar) {
    checkArgs(0, col);

    for (int i = 0; i < ROWS; ++i) {
        data[i][col] *= scalar;
    }
}

void Matrix::addMultipleRow(int targetRow, int sourceRow, double scalar) {
    checkArgs(targetRow, 0);
    checkArgs(sourceRow, 0);

    for (int i = 0; i < COLS; ++i) {
        data[targetRow][i] += data[sourceRow][i] * scalar;
    }
}

void Matrix::addMultipleCol(int targetCol, int sourceCol, double scalar) {
    checkArgs(0, targetCol);
    checkArgs(0, sourceCol);

    for (int i = 0; i < ROWS; ++i) {
        data[i][targetCol] += data[i][sourceCol] * scalar;
    }
}

int Matrix::getRows() const {
    return ROWS;
}

int Matrix::getCols() const {
    return COLS;
}

int *Matrix::getRow(int row) const {
    checkArgs(row, 0);

    int *result = new int[COLS];
    for (int i = 0; i < COLS; ++i) {
        result[i] = data[row][i];
    }
    return result;
}

int *Matrix::getCol(int col) const {
    checkArgs(0, col);

    int *result = new int[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        result[i] = data[i][col];
    }
    return result;
}

bool Matrix::isZeroRow(int row) const {
    checkArgs(row, 0);

    for (int i = 0; i < COLS; ++i) {
        if (data[row][i] != 0.0) {
            return false;
        }
    }
    return true;
}

bool Matrix::isZeroCol(int col) const {
    checkArgs(0, col);

    for (int i = 0; i < ROWS; ++i) {
        if (data[i][col] != 0.0) {
            return false;
        }
    }
    return true;
}

int Matrix::rank() const {
    Matrix rrefMatrix = rref();
    int rank = 0;

    while (rank <= ROWS) {
        if (rrefMatrix.isZeroRow(rank)) {
            return rank;
        }
        ++rank;
    };
    return rank;
}

const Matrix &Matrix::operator=(const Matrix &other) {
    if (&other == this) {
        return *this;
    }

    if (ROWS != other.ROWS || COLS != COLS) {
        throw std::invalid_argument("Matrix dimensions must match for assignment.");
    }

    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = other[i][j];
        }
    }
    return *this;
}

const Matrix &Matrix::operator=(const double **array) {
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = array[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator+(const Matrix &other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    }

    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] + other[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction.");
    }

    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] - other[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
    if (COLS != other.ROWS) {
        throw std::invalid_argument("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }

    Matrix result(ROWS, other.COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < other.COLS; ++j) {
            for (int k = 0; k < COLS; ++k) {
                result[i][j] += data[i][k] * other[k][j];
            }
        }
    }
    return result;
}

const Matrix &Matrix::operator+=(const Matrix &other) {
    return (*this = this->operator+(other));
}

const Matrix &Matrix::operator-=(const Matrix &other) {
    return (*this = this->operator-(other));
}

bool Matrix::operator==(const Matrix &other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        return false;
    }
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            if (data[i][j] != other[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix &other) const {
    return !(this->operator==(other));
}

bool Matrix::operator!() const {
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            if (data[i][j] != 0.0) {
                return false;
            }
        }
    }
    return true;
}

Matrix Matrix::transpose() const {
    Matrix result(COLS, ROWS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[j][i] = data[i][j];
        }
    }
    return result;
}

Matrix Matrix::ref() const {
    Matrix result(*this);
    int r = 0, r0 = 0, c = 0;
    double pivot;
    // double tmp;

    while (c < COLS) {
        while (r < ROWS && result[r][c] == 0.0) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }
        result.swapRows(r, r0);
        pivot = result[r0][c];
        for (int i = r0 + 1; i < ROWS; ++i) {
            if (result[i][c] != 0.0) {
                // tmp = result[i][c];
                // result. multiplyRow(i, pivot);
                // result.addMultipleRow(i, r0, -tmp);
                result.addMultipleRow(i, r0, -result[i][c] / pivot);
            }
        }
        r = ++r0;
        ++c;
    }
    return result;
}

Matrix Matrix::rref() const {
    Matrix result(*this);
    int r = 0, r0 = 0, c = 0;

    while (c < COLS) {
        while (r < ROWS && result[r][c] == 0.0) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }
        result.swapRows(r, r0);
        result.multiplyRow(r0, 1.0 / result[r0][c]);
        for (int i = 0; i < ROWS; ++i) {
            if (result[i][c] != 0.0 && i != r0) {
                result.addMultipleRow(i, r0, -result[i][c]);
            }
        }
        r = ++r0;
        ++c;
    }

    return result;
}

Matrix Matrix::operator()(int row, int col) const {
    checkArgs(row, col);

    Matrix result(ROWS - 1, COLS - 1);
    int i = 0, j = 0;
    for (int i = 0; i < ROWS - 1; ++i) {
        for (int j = 0; j < COLS - 1; ++j) {
            result[i][j] = data[i + (i >= row)][j + (j >= col)];
        }
    }
    return result;
}

double *Matrix::operator[](int row) {
    checkArgs(row, 0);
    return data[row];
}

const double *Matrix::operator[](int row) const {
    checkArgs(row, 0);
    return data[row];
}

Matrix::operator double**() const {
    return data;
}

Matrix operator*(double scalar, const Matrix &matrix) {
    Matrix result(matrix);
    for (int i = 0; i < result.ROWS; ++i) {
        result.multiplyCol(i, scalar);
    }
    return result;
}

std::ostream& operator<<(std::ostream &os, const Matrix &matrix) {
    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = 0; j < matrix.COLS; ++j) {
            os << matrix[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}

std::istream& operator>>(std::istream &is, Matrix &matrix) {
    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = 0; j < matrix.getCols(); ++j) {
            is >> matrix[i][j];
        }
    }
    return is;
}
