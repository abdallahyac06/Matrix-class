#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "matrix.hpp"
#include <iostream>

using std::ostream;
using std::istream;

template <typename T>
class SquareMatrix : public Matrix<T> {
    private:
        T determinantRecursive();

    public:
        SquareMatrix(unsigned long size = 1, const T& zero = T());
        SquareMatrix(const SquareMatrix<T>& other);
        SquareMatrix(SquareMatrix<T>&& other);
        SquareMatrix(const Matrix<T>& other);
        SquareMatrix(Matrix<T>&& other);
        SquareMatrix(const T* const* values, unsigned long size, const T& zero = T());
        virtual ~SquareMatrix() = default;

        SquareMatrix<T> transpose() const;
        SquareMatrix<T> ref() const;
        SquareMatrix<T> rref() const;
        SquareMatrix<T>& operator=(const SquareMatrix<T>& other);
        SquareMatrix<T> operator+(const Matrix<T>& other) const;
        SquareMatrix<T> operator-(const Matrix<T>& other) const;
        SquareMatrix<T> operator-() const;
        SquareMatrix<T> operator*(const SquareMatrix<T>& other) const;
        Matrix<T> operator*(const Matrix<T>& other) const;
        SquareMatrix<T> operator*(T scalar) const;
        SquareMatrix<T> operator/(T scalar) const;
        SquareMatrix<T>& operator+=(const Matrix<T>& other);
        SquareMatrix<T>& operator-=(const Matrix<T>& other);
        SquareMatrix<T>& operator*=(const Matrix<T>& other);
        SquareMatrix<T>& operator*=(T scalar);
        SquareMatrix<T>& operator/=(T scalar);

        static const SquareMatrix<T> id(unsigned long size, T scalar = T(1));
        T trace() const;
        T determinant() const;
        bool isLTriangular() const;
        bool isUTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;
        SquareMatrix<T> operator()(unsigned long row, unsigned long col) const;
        SquareMatrix<T> adjoint() const;
        SquareMatrix<T> inverse() const;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const SquareMatrix<T>& matrix);
template <typename T>
std::istream& operator>>(std::istream& is, SquareMatrix<T>& matrix);
template <typename T>
std::ostream& operator<<(std::ostream& os, const SquareMatrix<T>& matrix);

#include "squarematrix.tpp"

#endif