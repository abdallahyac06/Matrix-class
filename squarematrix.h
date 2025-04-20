#ifndef SQUAREMATRIX
#define SQUAREMATRIX

#include "matrix.h"
#include <iostream>

/**
 * @class SquareMatrix
 * @brief A class representing a square matrix, derived from the Matrix class.
 *
 * This class provides additional functionality specific to square matrices,
 * such as calculating the trace, determinant, and checking for special matrix
 * properties (e.g., triangular, diagonal, symmetric). It also supports matrix
 * inversion and adjoint computation.
 *
 * @note The class assumes that the input matrices are valid square matrices.
 */
class SquareMatrix : public Matrix {
    public:
        /**
         * @brief Constructs a SquareMatrix with a given size.
         * @param size The size of the square matrix (default is 1).
         */
        SquareMatrix(int size = 1);

        /**
         * @brief Copy constructor for SquareMatrix.
         * @param other The SquareMatrix to copy from.
         */
        SquareMatrix(const SquareMatrix &other);

        /**
         * @brief Move constructor for SquareMatrix.
         * @param other The SquareMatrix to move from.
         */
        SquareMatrix(SquareMatrix &&other);

        /**
         * @brief Constructs a SquareMatrix from a general Matrix.
         * @param matrix The Matrix to convert to a SquareMatrix.
         * @throws std::invalid_argument if the input matrix is not square.
         */
        SquareMatrix(const Matrix &matrix);

        /**
         * @brief Move constructor to create a SquareMatrix from a general Matrix.
         * @param matrix The Matrix to move and convert to a SquareMatrix.
         * @throws std::invalid_argument if the input matrix is not square.
         */
        SquareMatrix(Matrix &&matrix);

        /**
         * @brief Constructs a SquareMatrix from a 2D array.
         * @param array A pointer to a 2D array representing the matrix.
         * @param size The size of the square matrix.
         */
        SquareMatrix(double **array, int size);

        /**
         * @brief Calculates the trace of the square matrix.
         * @return The sum of the diagonal elements of the matrix.
         */
        double trace() const;

        /**
         * @brief Calculates the determinant of the square matrix.
         * @return The determinant of the matrix.
         */
        double determinant() const;

        /**
         * @brief Calculates the determinant of the square matrix using a recursive method.
         * @return The determinant of the matrix.
         */
        double determinantRecursive() const;

        /**
         * @brief Checks if the square matrix is lower triangular.
         * @return True if the matrix is lower triangular, false otherwise.
         */
        bool isLTriangular() const;

        /**
         * @brief Checks if the square matrix is upper triangular.
         * @return True if the matrix is upper triangular, false otherwise.
         */
        bool isUTriangular() const;

        /**
         * @brief Checks if the square matrix is diagonal.
         * @return True if the matrix is diagonal, false otherwise.
         */
        bool isDiagonal() const;

        /**
         * @brief Checks if the square matrix is symmetric.
         * @return True if the matrix is symmetric, false otherwise.
         */
        bool isSymmetric() const;

        /**
         * @brief Computes the adjoint of the square matrix.
         * @return A new SquareMatrix representing the adjoint of the matrix.
         */
        SquareMatrix adjoint() const;

        /**
         * @brief Computes the inverse of the square matrix.
         * @return A new SquareMatrix representing the inverse of the matrix.
         * @throws std::logic_error if the matrix is singular (non-invertible).
         */
        SquareMatrix inverse() const;

    /**
    * @brief Overloads the multiplication operator to scale a SquareMatrix by a scalar.
    * @param scalar The scalar value to multiply with.
    * @param matrix The SquareMatrix to scale.
    * @return A new SquareMatrix representing the scaled matrix.
    */
    friend SquareMatrix operator*(double scalar, const SquareMatrix &matrix);

    /**
     * @brief Overloads the stream insertion operator for SquareMatrix.
     * @param os The output stream.
     * @param matrix The SquareMatrix to output.
     * @return The output stream with the matrix data.
     */
    friend std::ostream &operator<<(std::ostream &os, const SquareMatrix &matrix);

    /**
     * @brief Overloads the stream extraction operator for SquareMatrix.
     * @param is The input stream.
     * @param matrix The SquareMatrix to populate with input data.
     * @return The input stream after reading the matrix data.
     */
    friend std::istream &operator>>(std::istream &is, SquareMatrix &matrix);
};

#endif