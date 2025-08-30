#include "matrix.h"
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>

using std::string;
using std::runtime_error;
using std::stringstream;
using std::ostream;
using std::istream;

using std::move;
using std::max;
using std::swap;
using std::setw;
using std::endl;

MatrixException::MatrixException(const string& message): runtime_error(message) {}

Matrix::Matrix(unsigned long rows, unsigned long cols): ROWS(rows), COLS(cols), data(ROWS, nullptr) {
    if (ROWS < 1 || COLS < 1) {
        throw MatrixException("Matrix dimensions must be greater than 0.");
    }
    
    for (double*& row: data) {
        row = new double[COLS];
        for (unsigned long i = 0; i < COLS; ++i) {
            row[i] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix& other): Matrix(other.data.data(), other.ROWS, other.COLS) {}

Matrix::Matrix(Matrix&& other): ROWS(other.ROWS), COLS(other.COLS) {
    data = move(other.data);
    other.data.assign(ROWS, nullptr);
}

Matrix::Matrix(const double* const* values, unsigned long rows, unsigned long cols): Matrix(rows, cols) {
    for (unsigned long i = 0; i < ROWS; ++i) {
        setRow(i, values[i]);
    }
}

Matrix::~Matrix() {
    for (double*& row: data) {
        delete[] row;
        row = nullptr;
    }
}

size_t Matrix::maxLength() const {
    size_t maxl = 0;
    stringstream ss;
    for (const double* row: data) {
        for (unsigned long i = 0; i < COLS; ++i) {
            ss.str("");
            ss << row[i];
            maxl = max(maxl, ss.str().length());
        }
    }
    
    return maxl;
}

void Matrix::setRow(unsigned long row, const double* values) {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    for (unsigned long i = 0; i < COLS; ++i) {
        data[row][i] = values[i];
    }
}

void Matrix::multiplyRow(unsigned long row, double scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[row][i] *= scalar;
    }
}

void Matrix::divideRow(unsigned long row, double scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[row][i] /= scalar;
    }
}

void Matrix::addMultipleRow(unsigned long targetRow, unsigned long sourceRow, double scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[targetRow][i] += data[sourceRow][i] * scalar;
    }
}

unsigned long Matrix::getRows() const {
    return ROWS;
}

unsigned long Matrix::getCols() const {
    return COLS;
}

bool Matrix::isZeroRow(unsigned long row) const {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    for (unsigned long i = 0; i < COLS; ++i) {
        if (data[row][i]) {
            return false;
        }
    }
    
    return true;
}

bool Matrix::isZeroCol(unsigned long col) const {
    if (col < 0 || col >= COLS) {
        throw MatrixException("Column index out of bounds.");
    }

    for (const double* row: data) {
        if (row[col]) {
            return false;
        }
    }
    
    return true;
}

unsigned long Matrix::rank() const {
    Matrix refMatrix(ref());
    unsigned long rank = 0;

    while (rank < ROWS && !refMatrix.isZeroRow(rank)) {
        ++rank;
    }

    return rank;
}

Matrix Matrix::transpose() const {
    Matrix result(COLS, ROWS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[j][i] = data[i][j];
        }
    }
    
    return result;
}

Matrix Matrix::ref() const {
    Matrix result(*this);
    unsigned long r = 0, r0 = 0, c = 0;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && !result[r][c]) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }

        swap(result.data[r], result.data[r0]);
        for (unsigned long i = r0 + 1; i < ROWS; ++i) {
            if (result[i][c]) {
                result.addMultipleRow(i, r0, -result[i][c] / result[r0][c], c + 1);
                result[i][c] = 0.0;
            }
        }
        r = ++r0;
        ++c;
    }
    
    return result;
}

Matrix Matrix::rref() const {
    Matrix result(*this);
    unsigned long r = 0, r0 = 0, c = 0;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && !result[r][c]) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }

        swap(result.data[r], result.data[r0]);
        result.divideRow(r0, result[r0][c], c);
        for (unsigned long i = 0; i < ROWS; ++i) {
            if (result[i][c] && i != r0) {
                result.addMultipleRow(i, r0, -result[i][c], c);
            }
        }
        r = ++r0;
        ++c;
    }

    return result;
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (&other == this) {
        return *this;
    }

    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw MatrixException("Matrix dimensions must match for assignment.");
    }

    for (unsigned long i = 0; i < ROWS; ++i) {
        setRow(i, other[i]);
    }
    
    return *this;
}

Matrix Matrix::operator+(const Matrix&other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw MatrixException("Matrix dimensions must match for addition.");
    }

    Matrix result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] + other[i][j];
        }
    }
    
    return result;
}

Matrix Matrix::operator-(const Matrix&other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw MatrixException("Matrix dimensions must match for subtraction.");
    }

    Matrix result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] - other[i][j];
        }
    }
    
    return result;
}

Matrix Matrix::operator-() const {
    Matrix result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = -data[i][j];
        }
    }
    
    return result;
}

Matrix Matrix::operator*(const Matrix&other) const {
    if (COLS != other.ROWS) {
        throw MatrixException("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }

    Matrix result(ROWS, other.COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < other.COLS; ++j) {
            result[i][j] = double();
            for (unsigned long k = 0; k < COLS; ++k) {
                result[i][j] += data[i][k] * other[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(*this);
    for (unsigned long i = 0; i < ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }
    
    return result;
}

Matrix Matrix::operator/(double scalar) const {
    Matrix result(*this);
    for (unsigned long i = 0; i < ROWS; ++i) {
        result.divideRow(i, scalar);
    }
    
    return result;
}

Matrix& Matrix::Matrix::operator+=(const Matrix&other) {
    return operator=(operator+(other));
}

Matrix& Matrix::operator-=(const Matrix&other) {
    return operator=(operator-(other));
}

Matrix& Matrix::operator*=(const Matrix&other) {
    return operator=(operator*(other));
}

Matrix& Matrix::operator*=(double scalar) {
    return operator=(operator*(scalar));
}

Matrix& Matrix::operator/=(double scalar) {
    return operator=(operator/(scalar));
}

bool Matrix::operator==(const Matrix&other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        return false;
    }

    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            if (data[i][j] != other[i][j]) {
                return false;
            }
        }
    }
    
    return true;
}

bool Matrix::operator!=(const Matrix&other) const {
    return !(operator==(other));
}

bool Matrix::operator!() const {
    return !(operator bool());
}

double*& Matrix::operator[](unsigned long row) {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }

    return data[row];
}

const double* Matrix::operator[](unsigned long row) const {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }

    return data[row];
}

Matrix::operator bool() const {
    for (unsigned long i = 0; i < ROWS; ++i) {
        if (!isZeroRow(i)) {
            return true;
        }
    }
    
    return false;
}

ostream& operator<<(ostream& os, const Matrix& matrix) {
    unsigned long maxl = 1 + matrix.maxLength();
    for (const double* row: matrix.data) {
        for (unsigned long i = 0; i < matrix.COLS; ++i) {
            os << setw(maxl) << row[i];
        }
        os << endl;
    }

    return os; 
}

istream& operator>>(istream& is, Matrix& matrix) {
    for (double* row: matrix.data) {
        for (unsigned long i = 0; i < matrix.COLS; ++i) {
            is >> row[i];
        }
    }
    
    return is;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}