#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>

using std::vector;
using std::ostream;
using std::istream;
using std::runtime_error;
using std::string;

template <typename T>
class Matrix {
    private:
        const int ROWS;
        const int COLS;
        
    protected:
        const T ZERO;
        vector<T*> data;
        size_t maxLength() const;
        void setRow(int row, const T* values);
        void multiplyRow(int row, T scalar, int startCol = 0);
        void divideRow(int row, T scalar, int startCol = 0);
        void addMultipleRow(int targetRow, int sourceRow, T scalar, int startCol = 0);

    public:
        Matrix(int rows = 1, int cols = 1, const T& zero = T());
        Matrix(const Matrix<T>& other);
        Matrix(Matrix<T>&& other);
        Matrix(const T* const* values, int rows = 1, int cols = 1, const T& zero = T());
        ~Matrix();

        int getRows() const;
        int getCols() const;
        void setRow(int row, const vector<T>& values);
        void setCol(int col, const vector<T>& values);
        vector<T> getRow(int row) const;
        vector<T> getCol(int col) const;
        bool isZeroRow(int row) const;
        bool isZeroCol(int col) const;
        int rank() const;
        Matrix<T> transpose() const;
        Matrix<T> ref() const;
        Matrix<T> rref() const;

        Matrix<T>& operator=(const Matrix<T>& other);
        Matrix<T> operator+(const Matrix<T>& other) const;
        Matrix<T> operator-(const Matrix<T>& other) const;
        Matrix<T> operator-() const;
        Matrix<T> operator*(const Matrix<T>& other) const;
        Matrix<T> operator*(T scalar) const;
        Matrix<T> operator/(T scalar) const;
        Matrix<T>& operator+=(const Matrix<T>& other);
        Matrix<T>& operator-=(const Matrix<T>& other);
        Matrix<T>& operator*=(const Matrix<T>& other);
        Matrix<T>& operator*=(T scalar);
        Matrix<T>& operator/=(T scalar);

        bool operator==(const Matrix<T>& other) const;
        bool operator!=(const Matrix<T>& other) const;
        bool operator!() const;
        T*& operator[](int row);
        const T* operator[](int row) const;
        operator bool() const;
    
    template <typename U>
    friend ostream& operator<<(ostream& os, const Matrix<U>& matrix); 
    template <typename U>
    friend istream& operator>>(istream& is, Matrix<U>& matrix);
};

class MatrixException : public runtime_error {
    public:
        explicit MatrixException(const string& message);
};

#include "matrix.tpp"

#endif