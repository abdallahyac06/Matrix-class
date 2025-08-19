#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <stdexcept>

using std::vector;
using std::ostream;
using std::istream;
using std::runtime_error;
using std::string;

class Matrix {
    private:
        const int ROWS;
        const int COLS;
        
    protected:
        vector<double*> data;
        size_t maxLength() const;
        void setRow(int row, const double* values);
        void multiplyRow(int row, double scalar, int startCol = 0);
        void divideRow(int row, double scalar, int startCol = 0);
        void addMultipleRow(int targetRow, int sourceRow, double scalar, int startCol = 0);

    public:
        Matrix(int rows = 1, int cols = 1);
        Matrix(const Matrix& other);
        Matrix(Matrix&& other);
        Matrix(const double* const* values, int rows = 1, int cols = 1);
        ~Matrix();

        int getRows() const;
        int getCols() const;
        void setRow(int row, const vector<double>& values);
        void setCol(int col, const vector<double>& values);
        vector<double> getRow(int row) const;
        vector<double> getCol(int col) const;
        bool isZeroRow(int row) const;
        bool isZeroCol(int col) const;
        int rank() const;
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
        double*& operator[](int row);
        const double* operator[](int row) const;
        operator bool() const;
        
    friend Matrix operator*(double scalar, const Matrix& matrix);
    friend ostream& operator<<(ostream& os, const Matrix& matrix); 
    friend istream& operator>>(istream& is, Matrix& matrix);
};

class MatrixException : public runtime_error {
    public:
        explicit MatrixException(const string& message);
};

#endif