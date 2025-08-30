#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::ostream;
using std::istream;
using std::runtime_error;
using std::string;

template <typename T>
class Matrix {
    protected:
        const unsigned long ROWS;
        const unsigned long COLS;
        const T ZERO;
        vector<T*> data;
        size_t maxLength() const;
        void setRow(unsigned long row, const T* values);
        void multiplyRow(unsigned long row, T scalar, unsigned long startCol = 0);
        void divideRow(unsigned long row, T scalar, unsigned long startCol = 0);
        void addMultipleRow(unsigned long targetRow, unsigned long sourceRow, T scalar, unsigned long startCol = 0);

    public:
        Matrix(unsigned long rows = 1, unsigned long cols = 1, const T& zero = T());
        Matrix(const Matrix<T>& other);
        Matrix(Matrix<T>&& other);
        Matrix(const T* const* values, unsigned long rows = 1, unsigned long cols = 1, const T& zero = T());
        ~Matrix();

        unsigned long getRows() const;
        unsigned long getCols() const;
        bool isZeroRow(unsigned long row) const;
        bool isZeroCol(unsigned long col) const;
        unsigned long rank() const;
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
        T*& operator[](unsigned long row);
        const T* operator[](unsigned long row) const;
        operator bool() const;
    
    template <typename U>
    friend ostream& operator<<(ostream& os, const Matrix<U>& matrix); 
    template <typename U>
    friend istream& operator>>(istream& is, Matrix<U>& matrix);
};

template <typename T>
ostream& operator<<(ostream& os, const Matrix<T>& matrix);

class MatrixException : public runtime_error {
    public:
        explicit MatrixException(const string& message);
};

#include "matrix.tpp"

#endif