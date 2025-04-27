#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <sstream>

using std::vector;
using std::ostream;
using std::istream;

template<typename T, T ZERO = T()>
class Matrix {
    private:
        const int ROWS;
        const int COLS;
        
    protected:
        vector<T*> data;

        size_t maxLength() const {
            size_t maxl = 0;
            std::stringstream ss;
            for (const T *row: data) {
                for (int i = 0; i < COLS; ++i) {
                    ss.str("");
                    ss << row[i];
                    maxl = std::max(maxl, ss.str().length());
                }
            }
            
            return maxl;
        }

        void setRow(int row, const T *values) {
            if (row < 0 || row >= ROWS) {
                throw std::out_of_range("Row index out of bounds.");
            }
            
            for (int i = 0; i < COLS; ++i) {
                data[row][i] = values[i];
            }
        }
        
        void multiplyRow(int row, T scalar, int startCol = 0) {
            for (int i = startCol; i < COLS; ++i) {
                data[row][i] *= scalar;
            }
        }
        
        void divideRow(int row, T scalar, int startCol = 0) {
            for (int i = startCol; i < COLS; ++i) {
                data[row][i] /= scalar;
            }
        }
        
        void addMultipleRow(int targetRow, int sourceRow, T scalar, int startCol = 0) {
            for (int i = startCol; i < COLS; ++i) {
                data[targetRow][i] += data[sourceRow][i] * scalar;
            }
        }

    public:
        Matrix(int rows, int cols): ROWS(rows), COLS(cols), data(ROWS, nullptr) {
            if (ROWS < 1 || COLS < 1) {
                throw std::invalid_argument("Matrix dimensions must be greater than 0.");
            }
            
            for (T *&row: data) {
                row = new T[COLS];
            }
        }

        Matrix(const Matrix<T, ZERO> &other): Matrix(other.data, other.ROWS, other.COLS) {}

        Matrix(Matrix<T, ZERO> &&other): ROWS(other.ROWS), COLS(other.COLS) {
            data = std::move(other.data);
            other.data.assign(ROWS, nullptr);
        }

        Matrix(const vector<T*> &values, int rows, int cols): Matrix(rows, cols) {
            for (int i = 0; i < ROWS; ++i) {
                setRow(i, values[i]);
            }
        }

        ~Matrix() {
            for (T *&ptr: data) {
                delete[] ptr;
            }
        }


        int getRows() const {
            return ROWS;
        }

        int getCols() const {
            return COLS;
        }

        void setRow(int row, const vector<T> &values) {
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

        void setCol(int col, const vector<T> &values) {
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

        vector<T> getRow(int row) const {
            if (row < 0 || row >= ROWS) {
                throw std::out_of_range("Row index out of bounds.");
            }
            
            return vector<T>(data[row], data[row] + COLS);
        }

        vector<T> getCol(int col) const {
            if (col < 0 || col >= COLS) {
                throw std::out_of_range("Column index out of bounds.");
            }
            
            vector<T> result(ROWS);
            for (int i = 0; i < ROWS; ++i) {
                result[i] = data[i][col];
            }
            
            return result;
        }

        bool isZeroRow(int row) const {
            if (row < 0 || row >= ROWS) {
                throw std::out_of_range("Row index out of bounds.");
            }
            
            for (int i = 0; i < COLS; ++i) {
                if (data[row][i] != ZERO) {
                    return false;
                }
            }
            
            return true;
        }

        bool isZeroCol(int col) const {
            if (col < 0 || col >= COLS) {
                throw std::out_of_range("Column index out of bounds.");
            }

            for (const T *row: data) {
                if (row[col] != ZERO) {
                    return false;
                }
            }
            
            return true;
        }

        int rank() const {
            Matrix<T, ZERO> refMatrix(ref());
            int rank = 0;

            while (rank < ROWS && !refMatrix.isZeroRow(rank)) {
                ++rank;
            }

            return rank;
        }

        const Matrix<T, ZERO> transpose() const {
            Matrix<T, ZERO> result(COLS, ROWS);
            for (int i = 0; i < ROWS; ++i) {
                for (int j = 0; j < COLS; ++j) {
                    result[j][i] = data[i][j];
                }
            }
            
            return result;
        }

        const Matrix<T, ZERO> ref() const {
            Matrix<T, ZERO> result(*this);
            int r = 0, r0 = 0, c = 0;
            T pivot;

            while (r0 < ROWS && c < COLS) {
                while (r < ROWS && result[r][c] == ZERO) {
                    ++r;
                }
                if (r == ROWS) {
                    r = r0;
                    ++c;
                    continue;
                }

                std::swap(result[r], result[r0]);
                for (int i = r0 + 1; i < ROWS; ++i) {
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

        const Matrix<T, ZERO> rref() const {
            Matrix<T, ZERO> result(*this);
            int r = 0, r0 = 0, c = 0;

            while (r0 < ROWS && c < COLS) {
                while (r < ROWS && result[r][c] == ZERO) {
                    ++r;
                }
                if (r == ROWS) {
                    r = r0;
                    ++c;
                    continue;
                }

                std::swap(result[r], result[r0]);
                result.divideRow(r0, result[r0][c], c);
                for (int i = 0; i < ROWS; ++i) {
                    if (result[i][c] == ZERO && i != r0) {
                        result.addMultipleRow(i, r0, -result[i][c]);
                    }
                }
                r = ++r0;
                ++c;
            }

            return result;
        }

        Matrix<T, ZERO> &operator=(const Matrix<T, ZERO> &other) {
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

        const Matrix<T, ZERO> operator+(const Matrix<T, ZERO> &other) const {
            if (ROWS != other.ROWS || COLS != other.COLS) {
                throw std::logic_error("Matrix dimensions must match for addition.");
            }

            Matrix<T, ZERO> result(ROWS, COLS);
            for (int i = 0; i < ROWS; ++i) {
                for (int j = 0; j < COLS; ++j) {
                    result[i][j] = data[i][j] + other[i][j];
                }
            }
            
            return result;
        }

        const Matrix<T, ZERO> operator-(const Matrix<T, ZERO> &other) const {
            if (ROWS != other.ROWS || COLS != other.COLS) {
                throw std::logic_error("Matrix dimensions must match for subtraction.");
            }

            Matrix<T, ZERO> result(ROWS, COLS);
            for (int i = 0; i < ROWS; ++i) {
                for (int j = 0; j < COLS; ++j) {
                    result[i][j] = data[i][j] - other[i][j];
                }
            }
            
            return result;
        }

        const Matrix<T, ZERO> operator-() const {
            Matrix<T, ZERO> result(ROWS, COLS);
            for (int i = 0; i < ROWS; ++i) {
                for (int j = 0; j < COLS; ++j) {
                    result[i][j] = -data[i][j];
                }
            }
            
            return result;
        }

        const Matrix<T, ZERO> operator*(const Matrix<T, ZERO> &other) const {
            if (COLS != other.ROWS) {
                throw std::logic_error("Number of columns in the first matrix must match the number of rows in the second matrix for multiplication.");
            }

            Matrix<T, ZERO> result(ROWS, other.COLS);
            for (int i = 0; i < ROWS; ++i) {
                for (int j = 0; j < other.COLS; ++j) {
                    result[i][j] = ZERO;
                    for (int k = 0; k < COLS; ++k) {
                        result[i][j] += data[i][k] * other[k][j];
                    }
                }
            }

            return result;
        }

        const Matrix<T, ZERO> operator*(T scalar) const {
            Matrix<T, ZERO> result(*this);
            for (int i = 0; i < ROWS; ++i) {
                result.multiplyRow(i, scalar);
            }
            
            return result;
        }

        const Matrix<T, ZERO> operator/(T scalar) const {
            Matrix<T, ZERO> result(*this);
            for (int i = 0; i < ROWS; ++i) {
                result.divideRow(i, scalar);
            }
            
            return result;
        }

        Matrix<T, ZERO> &operator+=(const Matrix<T, ZERO> &other) {
            return this->operator=(this->operator+(other));
        }

        Matrix<T, ZERO> &operator-=(const Matrix<T, ZERO> &other) {
            return this->operator=(this->operator-(other));
        }

        Matrix<T, ZERO> &operator*=(const Matrix<T, ZERO> &other) {
            return this->operator=(this->operator*(other));
        }

        Matrix<T, ZERO> &operator*=(T scalar) {
            return this->operator=(this->operator*(scalar));
        }

        Matrix<T, ZERO> &operator/=(T scalar) {
            return this->operator=(this->operator/(scalar));
        }

        bool operator==(const Matrix<T, ZERO> &other) const {
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

        bool operator!=(const Matrix<T, ZERO> &other) const {
            return !(this->operator==(other));
        }

        bool operator!() const {
            return !(this->operator bool());
        }

        T *&operator[](int row) {
            if (row < 0 || row >= ROWS) {
                throw std::out_of_range("Row index out of bounds.");
            }
            
            return data[row];
        }

        const T *operator[](int row) const {
            if (row < 0 || row >= ROWS) {
                throw std::out_of_range("Row index out of bounds.");
            }
            
            return data[row];
        }

        operator bool() const {
            for (int i = 0; i < ROWS; ++i) {
                if (!isZeroRow(i)) {
                    return true;
                }
            }
            
            return false;
        }
        
    friend Matrix<T, ZERO> operator*(T scalar, const Matrix<T, ZERO> &matrix) {
        return matrix * scalar;
    }

    friend ostream& operator<<(ostream &os, const Matrix<T, ZERO> &matrix) {
        int maxl = 1 + matrix.maxLength();
        for (const T *row: matrix.data) {
            for (int i = 0; i < matrix.COLS; ++i) {
                os << std::setw(maxl) << row[i];
            }
            os << std::endl;
        }

        return os;
    }

    friend istream& operator>>(istream &is, Matrix<T, ZERO> &matrix) {
        for (T *row: matrix.data) {
            for (int i = 0; i < matrix.COLS; ++i) {
                is >> row[i];
            }
        }
        
        return is;
    }
};