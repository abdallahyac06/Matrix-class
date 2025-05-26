#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "matrix.hpp"
#include <iostream>
#include <vector>

template <typename T>
class SquareMatrix : public Matrix<T> {
    private:
        T determinantRecursive();

    public:
        SquareMatrix(int size = 4, const T &zero = T());
        SquareMatrix(const SquareMatrix<T> &other);
        SquareMatrix(SquareMatrix<T> &&other);
        SquareMatrix(const Matrix<T> &other);
        SquareMatrix(Matrix<T> &&other);
        SquareMatrix(vector<T*> values, int size, const T &zero = T());
        virtual ~SquareMatrix() = default;

        static const SquareMatrix<T> id(int size);
        T trace() const;
        T determinant() const;
        bool isLTriangular() const;
        bool isUTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        SquareMatrix<T> operator()(int row, int col) const;
        SquareMatrix<T> adjoint() const;
        SquareMatrix<T> inverse() const;
        SquareMatrix<T> &operator=(const SquareMatrix<T> &other);

    template <typename U>
    friend const SquareMatrix<U> operator*(U scalar, const SquareMatrix<U> &matrix);
    template <typename U>
    friend ostream &operator<<(ostream &os, const SquareMatrix<U> &matrix);
    template <typename U>
    friend istream &operator>>(istream &is, SquareMatrix<U> &matrix);
};

#include "squarematrix.tpp"

#endif