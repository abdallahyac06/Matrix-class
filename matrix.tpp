#include "matrix.hpp"
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

template <typename T>
Matrix<T>::Matrix(unsigned long rows, unsigned long cols, const T& zero): ROWS(rows), COLS(cols), ZERO(zero), data(ROWS, nullptr) {
    if (ROWS < 1 || COLS < 1) {
        throw MatrixException("Matrix dimensions must be greater than 0.");
    }
    
    for (T*& row: data) {
        row = new T[COLS];
        for (unsigned long i = 0; i < COLS; ++i) {
            row[i] = ZERO;
        }
    }
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other): Matrix(other.data.data(), other.ROWS, other.COLS, other.ZERO) {}

template <typename T>
Matrix<T>::Matrix(Matrix<T>&& other): ROWS(other.ROWS), COLS(other.COLS), ZERO(other.ZERO) {
    data = move(other.data);
    other.data.assign(ROWS, nullptr);
}

template <typename T>
Matrix<T>::Matrix(const T* const* values, unsigned long rows, unsigned long cols, const T& zero): Matrix(rows, cols, zero) {
    for (unsigned long i = 0; i < ROWS; ++i) {
        setRow(i, values[i]);
    }
}

template <typename T>
Matrix<T>::~Matrix() {
    for (T*& row: data) {
        delete[] row;
        row = nullptr;
    }
}

template <typename T>
size_t Matrix<T>::maxLength() const {
    size_t maxl = 0;
    stringstream ss;
    for (const T* row: data) {
        for (unsigned long i = 0; i < COLS; ++i) {
            ss.str("");
            ss << row[i];
            maxl = max(maxl, ss.str().length());
        }
    }
    
    return maxl;
}

template <typename T>
void Matrix<T>::setRow(unsigned long row, const T* values) {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    for (unsigned long i = 0; i < COLS; ++i) {
        data[row][i] = values[i];
    }
}

template <typename T>
void Matrix<T>::multiplyRow(unsigned long row, T scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[row][i] *= scalar;
    }
}

template <typename T>
void Matrix<T>::divideRow(unsigned long row, T scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[row][i] /= scalar;
    }
}

template <typename T>
void Matrix<T>::addMultipleRow(unsigned long targetRow, unsigned long sourceRow, T scalar, unsigned long startCol) {
    for (unsigned long i = startCol; i < COLS; ++i) {
        data[targetRow][i] += data[sourceRow][i] * scalar;
    }
}

template <typename T>
unsigned long Matrix<T>::getRows() const {
    return ROWS;
}

template <typename T>
unsigned long Matrix<T>::getCols() const {
    return COLS;
}

template <typename T>
bool Matrix<T>::isZeroRow(unsigned long row) const {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    for (unsigned long i = 0; i < COLS; ++i) {
        if (data[row][i] != ZERO) {
            return false;
        }
    }
    
    return true;
}

template <typename T>
bool Matrix<T>::isZeroCol(unsigned long col) const {
    if (col < 0 || col >= COLS) {
        throw MatrixException("Column index out of bounds.");
    }

    for (const T* row: data) {
        if (row[col] != ZERO) {
            return false;
        }
    }
    
    return true;
}

template <typename T>
unsigned long Matrix<T>::rank() const {
    Matrix<T> refMatrix(ref());
    unsigned long rank = 0;

    while (rank < ROWS && !refMatrix.isZeroRow(rank)) {
        ++rank;
    }

    return rank;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> result(COLS, ROWS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[j][i] = data[i][j];
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::ref() const {
    Matrix<T> result(*this);
    unsigned long r = 0, r0 = 0, c = 0;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && result[r][c] == ZERO) {
            ++r;
        }
        if (r >= ROWS) {
            r = r0;
            ++c;
            continue;
        }

        swap(result[r], result[r0]);
        for (unsigned long i = r0 + 1; i < ROWS; ++i) {
            if (result[i][c] != ZERO) {
                result.addMultipleRow(i, r0, -result[i][c] / result[r0][c], c + 1);
                result[i][c] = ZERO;
            }
        }
        r = ++r0;
        ++c;
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::rref() const {
    Matrix<T> result(*this);
    unsigned long r = 0, r0 = 0, c = 0;

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && result[r][c] == ZERO) {
            ++r;
        }
        if (r >= ROWS) {
            r = r0;
            ++c;
            continue;
        }

        swap(result[r], result[r0]);
        result.divideRow(r0, result[r0][c], c);
        for (unsigned long i = 0; i < ROWS; ++i) {
            if (result[i][c] != ZERO && i != r0) {
                result.addMultipleRow(i, r0, -result[i][c]);
            }
        }
        r = ++r0;
        ++c;
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
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

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw MatrixException("Matrix dimensions must match for addition.");
    }

    Matrix<T> result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] + other[i][j];
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (ROWS != other.ROWS || COLS != other.COLS) {
        throw MatrixException("Matrix dimensions must match for subtraction.");
    }

    Matrix<T> result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = data[i][j] - other[i][j];
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> result(ROWS, COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < COLS; ++j) {
            result[i][j] = -data[i][j];
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (COLS != other.ROWS) {
        throw MatrixException("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }

    Matrix<T> result(ROWS, other.COLS);
    for (unsigned long i = 0; i < ROWS; ++i) {
        for (unsigned long j = 0; j < other.COLS; ++j) {
            result[i][j] = ZERO;
            for (unsigned long k = 0; k < COLS; ++k) {
                result[i][j] += data[i][k] * other[k][j];
            }
        }
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(T scalar) const {
    Matrix<T> result(*this);
    for (unsigned long i = 0; i < ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(T scalar) const {
    Matrix<T> result(*this);
    for (unsigned long i = 0; i < ROWS; ++i) {
        result.divideRow(i, scalar);
    }
    
    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
    return operator=(operator+(other));
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
    return operator=(operator-(other));
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
    return operator=(operator*(other));
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(T scalar) {
    return operator=(operator*(scalar));
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(T scalar) {
    return operator=(operator/(scalar));
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& other) const {
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

template <typename T>
bool Matrix<T>::operator!=(const Matrix<T>& other) const {
    return !(operator==(other));
}

template <typename T>
bool Matrix<T>::operator!() const {
    return !(operator bool());
}

template <typename T>
T*& Matrix<T>::operator[](unsigned long row) {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    return data[row];
}

template <typename T>
const T* Matrix<T>::operator[](unsigned long row) const {
    if (row < 0 || row >= ROWS) {
        throw MatrixException("Row index out of bounds.");
    }
    
    return data[row];
}

template <typename T>
Matrix<T>::operator bool() const {
    for (unsigned long i = 0; i < ROWS; ++i) {
        if (!isZeroRow(i)) {
            return true;
        }
    }
    
    return false;
}

template <typename T>
ostream& operator<<(ostream& os, const Matrix<T>& matrix) {
    unsigned long maxl = 1 + matrix.maxLength();
    for (const T* row: matrix.data) {
        for (unsigned long i = 0; i < matrix.COLS; ++i) {
            os << setw(maxl) << row[i];
        }
        os << endl;
    }

    return os;
}

template <typename T>
istream& operator>>(istream& is, Matrix<T>& matrix) {
    for (T* row: matrix.data) {
        for (unsigned long i = 0; i < matrix.COLS; ++i) {
            is >> row[i];
        }
    }

    return is;
}

template <typename T>
Matrix<T> operator*(T scalar, const Matrix<T>& matrix) {
    return matrix * scalar;
}