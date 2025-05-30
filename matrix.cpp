#include "matrix.h"
#include <stdexcept>
#include <iomanip>
#include <sstream>

Matrix::Matrix(int rows, int cols): ROWS(rows), COLS(cols), data(ROWS, nullptr) {
    if (ROWS < 1 || COLS < 1) {
        throw std::invalid_argument("Matrix dimensions must be greater than 0.");
    }
    
    for (double *&row: data) {
        row = new double[COLS];
        for (int i = 0; i < COLS; ++i) {
            row[i] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix&other): Matrix(other.data, other.ROWS, other.COLS) {}

Matrix::Matrix(Matrix &&other): ROWS(other.ROWS), COLS(other.COLS) {
    data = std::move(other.data);
    other.data.assign(ROWS, nullptr);
}

Matrix::Matrix(const vector<double*> &values, int rows, int cols): Matrix(rows, cols) {
    for (int i = 0; i < ROWS; ++i) {
        setRow(i, values[i]);
    }
}

Matrix::~Matrix() {
    for (double *&row: data) {
        delete[] row;
    }
}

size_t Matrix::maxLength() const {
    size_t maxl = 0;
    std::stringstream ss;
    for (const double *row: data) {
        for (int i = 0; i < COLS; ++i) {
            ss.str("");
            ss << row[i];
            maxl = std::max(maxl, ss.str().length());
        }
    }
    
    return maxl;
}

void Matrix::setRow(int row, const double *values) {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }
    
    for (int i = 0; i < COLS; ++i) {
        data[row][i] = values[i];
    }
}

void Matrix::multiplyRow(int row, double scalar, int startCol) {
    for (int i = startCol; i < COLS; ++i) {
        data[row][i] *= scalar;
    }
}

void Matrix::divideRow(int row, double scalar, int startCol) {
    for (int i = startCol; i < COLS; ++i) {
        data[row][i] /= scalar;
    }
}

void Matrix::addMultipleRow(int targetRow, int sourceRow, double scalar, int startCol) {
    for (int i = startCol; i < COLS; ++i) {
        data[targetRow][i] += data[sourceRow][i] * scalar;
    }
}

int Matrix::getRows() const {
    return ROWS;
}

int Matrix::getCols() const {
    return COLS;
}

void Matrix::setRow(int row, const vector<double> &values) {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }

    if (values.size() != COLS) {
        throw std::invalid_argument("The vector size must be equal to the number columns.");
    }
    
    for (int i = 0; i < COLS; ++i) {
        data[row][i] = values[i];
    }
}

void Matrix::setCol(int col, const vector<double> &values) {
    if (col < 0 || col >= COLS) {
        throw std::out_of_range("Column index out of bounds.");
    }

    if (values.size() != ROWS) {
        throw std::invalid_argument("The vector size must be equal to the number rows.");
    }

    for (int i = 0; i < ROWS; ++i) {
        data[i][col] = values[i];
    }
}

vector<double> Matrix::getRow(int row) const {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }

    return vector<double>(data[row], data[row] + COLS);
}

vector<double> Matrix::getCol(int col) const {
    if (col < 0 || col >= COLS) {
        throw std::out_of_range("Column index out of bounds.");
    }
    
    vector<double> result(ROWS);
    for (int i = 0; i < ROWS; ++i) {
        result[i] = data[i][col];
    }
    
    return result;
}

bool Matrix::isZeroRow(int row) const {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }
    
    for (int i = 0; i < COLS; ++i) {
        if (data[row][i]) {
            return false;
        }
    }
    
    return true;
}

bool Matrix::isZeroCol(int col) const {
    if (col < 0 || col >= COLS) {
        throw std::out_of_range("Column index out of bounds.");
    }

    for (const double *row: data) {
        if (row[col]) {
            return false;
        }
    }
    
    return true;
}

int Matrix::rank() const {
    Matrix refMatrix(ref());
    int rank = 0;

    while (rank < ROWS && !refMatrix.isZeroRow(rank)) {
        ++rank;
    }

    return rank;
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

    while (r0 < ROWS && c < COLS) {
        while (r < ROWS && !result[r][c]) {
            ++r;
        }
        if (r == ROWS) {
            r = r0;
            ++c;
            continue;
        }

        std::swap(result.data[r], result.data[r0]);
        for (int i = r0 + 1; i < ROWS; ++i) {
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

        std::swap(result.data[r], result.data[r0]);
        result.divideRow(r0, result[r0][c], c);
        for (int i = 0; i < ROWS; ++i) {
            if (result[i][c] && i != r0) {
                result.addMultipleRow(i, r0, -result[i][c], c);
            }
        }
        r = ++r0;
        ++c;
    }

    return result;
}

Matrix &Matrix::operator=(const Matrix&other) {
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

Matrix Matrix::operator+(const Matrix&other) const {
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

Matrix Matrix::operator-(const Matrix&other) const {
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

Matrix Matrix::operator-() const {
    Matrix result(ROWS, COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            result[i][j] = -data[i][j];
        }
    }
    
    return result;
}

Matrix Matrix::operator*(const Matrix&other) const {
    if (COLS != other.ROWS) {
        throw std::logic_error("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
    }

    Matrix result(ROWS, other.COLS);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < other.COLS; ++j) {
            result[i][j] = double();
            for (int k = 0; k < COLS; ++k) {
                result[i][j] += data[i][k] * other[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(*this);
    for (int i = 0; i < ROWS; ++i) {
        result.multiplyRow(i, scalar);
    }
    
    return result;
}

Matrix Matrix::operator/(double scalar) const {
    Matrix result(*this);
    for (int i = 0; i < ROWS; ++i) {
        result.divideRow(i, scalar);
    }
    
    return result;
}

Matrix &Matrix::Matrix::operator+=(const Matrix&other) {
    return operator=(operator+(other));
}

Matrix &Matrix::operator-=(const Matrix&other) {
    return operator=(operator-(other));
}

Matrix &Matrix::operator*=(const Matrix&other) {
    return operator=(operator*(other));
}

Matrix &Matrix::operator*=(double scalar) {
    return operator=(operator*(scalar));
}

Matrix &Matrix::operator/=(double scalar) {
    return operator=(operator/(scalar));
}

bool Matrix::operator==(const Matrix&other) const {
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

bool Matrix::operator!=(const Matrix&other) const {
    return !(operator==(other));
}

bool Matrix::operator!() const {
    return !(operator bool());
}

double *&Matrix::operator[](int row) {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }

    return data[row];
}

const double *Matrix::operator[](int row) const {
    if (row < 0 || row >= ROWS) {
        throw std::out_of_range("Row index out of bounds.");
    }

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
    return matrix * scalar;
}

ostream& operator<<(ostream &os, const Matrix &matrix) {
    int maxl = 1 + matrix.maxLength();
    os << std::left;
    for (const double *row: matrix.data) {
        for (int i = 0; i < matrix.COLS; ++i) {
            os << std::setw(maxl) << row[i];
        }
        os << std::endl;
    }

    return os << std::right;
}

istream& operator>>(istream &is, Matrix &matrix) {
    for (double *&row: matrix.data) {
        for (int i = 0; i < matrix.COLS; ++i) {
            is >> row[i];
        }
    }
    
    return is;
}
