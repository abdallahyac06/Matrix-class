#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::ostream;
using std::istream;
using std::runtime_error;
using std::string;

class Matrix {
    protected:
        const unsigned long ROWS;
        const unsigned long COLS;
        vector<double*> data;
        size_t maxLength() const;
        void setRow(unsigned long row, const double* values);
        void multiplyRow(unsigned long row, double scalar, unsigned long startCol = 0);
        void divideRow(unsigned long row, double scalar, unsigned long startCol = 0);
        void addMultipleRow(unsigned long targetRow, unsigned long sourceRow, double scalar, unsigned long startCol = 0);

    public:
        Matrix(unsigned long rows = 1, unsigned long cols = 1);
        Matrix(const Matrix& other);
        Matrix(Matrix&& other);
        Matrix(const double* const* values, unsigned long rows = 1, unsigned long cols = 1);
        ~Matrix();

        unsigned long getRows() const;
        unsigned long getCols() const;
        bool isZeroRow(unsigned long row) const;
        bool isZeroCol(unsigned long col) const;
        unsigned long rank() const;
        Matrix transpose() const;
        Matrix ref() const;
        Matrix rref() const;

        Matrix& operator=(const Matrix&other);
        Matrix operator+(const Matrix&other) const;
        Matrix operator-(const Matrix&other) const;
        Matrix operator-() const;
        Matrix operator*(const Matrix&other) const;
        Matrix operator*(double scalar) const;
        Matrix operator/(double scalar) const;
        Matrix& operator+=(const Matrix&other);
        Matrix& operator-=(const Matrix&other);
        Matrix& operator*=(const Matrix&other);
        Matrix& operator*=(double scalar);
        Matrix& operator/=(double scalar);
        bool operator==(const Matrix&other) const;
        bool operator!=(const Matrix&other) const;
        bool operator!() const;
        double*& operator[](unsigned long row);
        const double* operator[](unsigned long row) const;
        operator bool() const;
        
    friend ostream& operator<<(ostream& os, const Matrix& matrix); 
    friend istream& operator>>(istream& is, Matrix& matrix);
};

Matrix operator*(double scalar, const Matrix& matrix);

class MatrixException : public runtime_error {
    public:
        explicit MatrixException(const string& message);
};

#endif