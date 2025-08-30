#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "matrix.h"
#include <iostream>

using std::ostream;
using std::istream;

class SquareMatrix : public Matrix {
    private:
        double determinantRecursive();
        
    public:
        SquareMatrix(unsigned long size = 1);
        SquareMatrix(const SquareMatrix& other);
        SquareMatrix(SquareMatrix&& other);
        SquareMatrix(const Matrix& other);
        SquareMatrix(Matrix&& other);
        SquareMatrix(const double* const* values, unsigned long size = 1);
        virtual ~SquareMatrix() = default;

        SquareMatrix transpose() const;
        SquareMatrix ref() const;
        SquareMatrix rref() const;
        SquareMatrix& operator=(const SquareMatrix& other);
        SquareMatrix operator+(const Matrix& other) const;
        SquareMatrix operator-(const Matrix& other) const;
        SquareMatrix operator-() const;
        SquareMatrix operator*(const SquareMatrix& other) const;
        Matrix operator*(const Matrix& other) const;
        SquareMatrix operator*(double scalar) const;
        SquareMatrix operator/(double scalar) const;
        SquareMatrix& operator+=(const Matrix& other);
        SquareMatrix& operator-=(const Matrix& other);
        SquareMatrix& operator*=(const Matrix& other);
        SquareMatrix& operator*=(double scalar);
        SquareMatrix& operator/=(double scalar);

        static SquareMatrix id(unsigned long size);
        double trace() const;
        double determinant() const;
        bool isLowerTriangular() const;
        bool isUpperTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        SquareMatrix operator()(unsigned long row, unsigned long col) const;
        SquareMatrix adjoint() const;
        SquareMatrix inverse() const;
};

ostream& operator<<(ostream& os, const SquareMatrix& matrix);
istream& operator>>(istream& is, SquareMatrix& matrix);
SquareMatrix operator*(double scalar, const SquareMatrix& matrix);

#endif