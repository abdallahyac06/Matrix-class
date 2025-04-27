#include "matrix.h"
#include <iostream>
#include <vector>

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

class SquareMatrix : public Matrix {
    public:
        SquareMatrix(int size = 1);
        SquareMatrix(const SquareMatrix &other);
        SquareMatrix(SquareMatrix &&other);
        SquareMatrix(const Matrix &other);
        SquareMatrix(Matrix &&other);
        SquareMatrix(vector<double*> values, int size);
        virtual ~SquareMatrix() = default;

        double trace() const;
        double determinant() const;
        bool isLowerTriangular() const;
        bool isUpperTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        const SquareMatrix operator()(int row, int col) const;
        const SquareMatrix adjoint() const;
        const SquareMatrix inverse() const;
        SquareMatrix &operator=(const SquareMatrix &other);

    friend const SquareMatrix operator*(double scalar, const SquareMatrix &matrix);
    friend ostream &operator<<(ostream &os, const SquareMatrix &matrix);
    friend istream &operator>>(istream &is, SquareMatrix &matrix);
};

#endif