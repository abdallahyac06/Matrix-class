#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "matrix.h"

class SquareMatrix : public Matrix {
    private:
        double determinantRecursive();
        
    public:
        SquareMatrix(int size = 1);
        SquareMatrix(const SquareMatrix& other);
        SquareMatrix(SquareMatrix&& other);
        SquareMatrix(const Matrix& other);
        SquareMatrix(Matrix&& other);
        SquareMatrix(const double* const* values, int size = 1);
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

        static SquareMatrix id(int size);
        double trace() const;
        double determinant() const;
        bool isLowerTriangular() const;
        bool isUpperTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        SquareMatrix operator()(int row, int col) const;
        SquareMatrix adjoint() const;
        SquareMatrix inverse() const;
};

#endif