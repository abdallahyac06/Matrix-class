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
    }
}

Matrix::Matrix(const Matrix &other): Matrix(other.data, other.ROWS, other.COLS) {}

Matrix::Matrix(Matrix &&other): ROWS(other.ROWS), COLS(other.COLS) {
    data = other.data;
    other.data = nullptr;
}

Matrix::Matrix(double **array, int rows, int cols): Matrix(rows, cols) {
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            data[i][j] = array[i][j];
        }
    }
}

Matrix::~Matrix() {
    if (data) {
        for (int i = 0; i < ROWS; ++i) {
            delete[] data[i];
            data[i] = nullptr;
        }

        delete[] data;
        data = nullptr;
    }
}

void Matrix::checkArgs(int row, int col) const {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }
    if (col < 0 || col >= COLS) {
        throw std::out_of_range("Col index out of bounds.");
    }
}

void Matrix::swapRows(int row1, int row2) {
    checkArgs(row1, 0);
    checkArgs(row2, 0);

    std::swap(data[row1], data[row2]);
}

void Matrix::multiplyRow(int row, double scalar) {
    checkArgs(row, 0);

    for (int i = 0; i < COLS; ++i) {
        data[row][i] *= scalar;
    }
}

void Matrix::addMultipleRow(int targetRow, int sourceRow, double scalar) {
    checkArgs(targetRow, 0);
    checkArgs(sourceRow, 0);

    for (int i = 0; i < COLS; ++i) {
        data[targetRow][i] += data[sourceRow][i] * scalar;
    }
}

double **Matrix::getData() const {
    return data;
}

int Matrix::getRows() const {
    return ROWS;
}

int Matrix::getCols() const {
    return COLS;
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
        if (data[row][i]) {
            return false;
        }
    }
    
    return true;
}

bool Matrix::isZeroCol(int col) const {
    checkArgs(0, col);

    for (int i = 0; i < ROWS; ++i) {
        if (data[i][col]) {
            return false;
        }
    }
    
    return true;
}

int Matrix::rank() const {
    Matrix rrefMatrix(std::move(rref()));
    int rank = 0;

    while (rank < ROWS && !rrefMatrix.isZeroRow(rank)) {
        ++rank;
    }

    return rank;
}    

const Matrix Matrix::transpose() const {
    Matrix result(COLS, ROWS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[j][i] = data[i][j];
        }
    }
    
    return result;
}

const Matrix Matrix::ref() const {
    Matrix result(*this);
    int r = 0, r0 = 0, c = 0;
    double pivot;
    // double tmp;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && !result[r][c]) {
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
            if (result[i][c]) {
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

const Matrix Matrix::rref(double unit) const {
    Matrix result(*this);
    int r = 0, r0 = 0, c = 0;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && !result[r][c]) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }

        result.swapRows(r, r0);
        result.multiplyRow(r0, unit / result[r0][c]);
        for (int i = 0; i < ROWS; ++i) {
            if (result[i][c] && i != r0) {
                result.addMultipleRow(i, r0, -result[i][c]);
            }
        }
        r = ++r0;
        ++c;
    }

    return result;
}

Matrix &Matrix::operator=(const Matrix &other) {
    if (&other == this) {
        return *this;
    }

    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw std::logic_error("Matrix dimensions must match for assignment.");
    }    

    for (int i = 0; i < ROWS; ++i) {
        setRow(i, other[i]);
    }    
    
    return *this;
}    

const Matrix Matrix::operator+(const Matrix &other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw std::logic_error("Matrix dimensions must match for addition.");
    }    

    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] + other[i][j];
        }    
    }    
    
    return result;
}    

const Matrix Matrix::operator-(const Matrix &other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw std::logic_error("Matrix dimensions must match for subtraction.");
    }    

    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] - other[i][j];
        }    
    }    
    
    return result;
}    

const Matrix Matrix::operator*(const Matrix &other) const {
    if (COLS != other.ROWS) {
        throw std::logic_error("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }    

    Matrix result(ROWS, other.COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < other.COLS; ++j) {
            result[i][j] = 0.0;
            for (int k = 0; k < COLS; ++k) {
                result[i][j] += data[i][k] * other[k][j];
            }    
        }    
    }    
    
    return result;
}    

const Matrix Matrix::operator*(double scalar) const {
    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }    
    
    return result;
}    

Matrix &Matrix::operator+=(const Matrix &other) {
    return this->operator=(this->operator+(other));
}    

Matrix &Matrix::operator-=(const Matrix &other) {
    return this->operator=(this->operator-(other));
}    

Matrix &Matrix::operator*=(const Matrix &other) {
    return this->operator=(this->operator*(other));
}    

Matrix &Matrix::operator*=(double scalar) {
    return this->operator=(this->operator*(scalar));
}

const Matrix Matrix::operator-() const {
    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = -data[i][j];
        }    
    }
    
    return result;
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
    return !(this->operator bool());
}

double *&Matrix::operator[](int row) {
    checkArgs(row, 0);
    return data[row];
}

const double *Matrix::operator[](int row) const {
    checkArgs(row, 0);
    return data[row];
}

Matrix::operator bool() const {
    for (int i = 0; i < ROWS; ++i) {
        if (!isZeroRow(i)) {
            return true;
        }
    }
    
    return false;
}

Matrix operator*(double scalar, const Matrix &matrix) {
    Matrix result(matrix);
    for (int i = 0; i < result.ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }
    
    return result;
}

std::ostream& operator<<(std::ostream &os, const Matrix &matrix) {
    for (int i = 0; i < matrix.ROWS; ++i) {
        for (int j = 0; j < matrix.COLS; ++j) {
            os << matrix[i][j] << " ";
        }
        os << std::endl;
    }
    
    return os;
}

std::istream& operator>>(std::istream &is, Matrix &matrix) {
    for (int i = 0; i < matrix.ROWS; ++i) {
        for (int j = 0; j < matrix.COLS; ++j) {
            is >> matrix[i][j];
        }
    }
    
    return is;
}
