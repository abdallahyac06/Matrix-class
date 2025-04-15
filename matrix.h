#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

/**
 * @class Matrix
 * @brief A class for representing and manipulating matrices.
 * 
 * This class provides a variety of methods for matrix operations, including
 * basic arithmetic, row and column manipulations, and transformations such as
 * transpose, row echelon form (REF), and reduced row echelon form (RREF).
 */
class Matrix {
    private:
        /**
         * @brief Number of rows in the matrix.
         */
        const int ROWS;

        /**
         * @brief Number of columns in the matrix.
         */
        const int COLS;

        /**
         * @brief Pointer to the 2D array storing matrix data.
         */
        double** data;

        /**
         * @brief Validates the row and column indices.
         * @param row Row index to check.
         * @param col Column index to check.
         * @throws std::out_of_range if indices are invalid.
         */
        void checkArgs(int row, int col) const;

        /**
         * @brief Sets the values of a specific row.
         * @param row Row index to set.
         * @param values Array of values to assign to the row.
         */
        void setRow(int row, const double* values);

        /**
         * @brief Sets the values of a specific column.
         * @param col Column index to set.
         * @param values Array of values to assign to the column.
         */
        void setCol(int col, const double* values);

        /**
         * @brief Swaps two rows in the matrix.
         * @param row1 Index of the first row.
         * @param row2 Index of the second row.
         */
        void swapRows(int row1, int row2);

        /**
         * @brief Swaps two columns in the matrix.
         * @param col1 Index of the first column.
         * @param col2 Index of the second column.
         */
        void swapCols(int col1, int col2);

        /**
         * @brief Multiplies a row by a scalar value.
         * @param row Row index to multiply.
         * @param scalar Scalar value to multiply the row by.
         */
        void multiplyRow(int row, double scalar);

        /**
         * @brief Multiplies a column by a scalar value.
         * @param col Column index to multiply.
         * @param scalar Scalar value to multiply the column by.
         */
        void multiplyCol(int col, double scalar);

        /**
         * @brief Adds a multiple of one row to another row.
         * @param targetRow Index of the row to modify.
         * @param sourceRow Index of the row to scale and add.
         * @param scalar Scalar value to multiply the source row by.
         */
        void addMultipleRow(int targetRow, int sourceRow, double scalar);

        /**
         * @brief Adds a multiple of one column to another column.
         * @param targetCol Index of the column to modify.
         * @param sourceCol Index of the column to scale and add.
         * @param scalar Scalar value to multiply the source column by.
         */
        void addMultipleCol(int targetCol, int sourceCol, double scalar);

    public:
        /**
         * @brief Constructs a matrix with the specified dimensions.
         * @param rows Number of rows (default is 1).
         * @param cols Number of columns (default is 1).
         */
        Matrix(int rows = 1, int cols = 1);

        /**
         * @brief Copy constructor.
         * @param other Matrix to copy from.
         */
        Matrix(const Matrix &other);

        /**
         * @brief Constructs a matrix from a 2D array.
         * @param array Pointer to a 2D array of values.
         * @param rows Number of rows in the array.
         * @param cols Number of columns in the array.
         */
        Matrix(const double** array, int rows, int cols);

        /**
         * @brief Destructor to free allocated memory.
         */
        ~Matrix();

        /**
         * @brief Gets the number of rows in the matrix.
         * @return Number of rows.
         */
        int getRows() const;

        /**
         * @brief Gets the number of columns in the matrix.
         * @return Number of columns.
         */
        int getCols() const;

        /**
         * @brief Retrieves a specific row as an array.
         * @param row Index of the row to retrieve.
         * @return Pointer to the row array.
         */
        int* getRow(int row) const;

        /**
         * @brief Retrieves a specific column as an array.
         * @param col Index of the column to retrieve.
         * @return Pointer to the column array.
         */
        int* getCol(int col) const;

        /**
         * @brief Checks if a specific row is a zero row.
         * @param row Index of the row to check.
         * @return True if the row is a zero row, false otherwise.
         */
        bool isZeroRow(int row) const;

        /**
         * @brief Checks if a specific column is a zero column.
         * @param col Index of the column to check.
         * @return True if the column is a zero column, false otherwise.
         */
        bool isZeroCol(int col) const;

        /**
         * @brief Computes the rank of the matrix.
         * @return Rank of the matrix.
         */
        int rank() const;

        /**
         * @brief Assignment operator to copy another matrix.
         * @param other Matrix to copy from.
         * @return Reference to the current matrix.
         */
        const Matrix &operator=(const Matrix &other);

        /**
         * @brief Assignment operator to copy from a 2D array.
         * @param array Pointer to a 2D array of values.
         * @return Reference to the current matrix.
         */
        const Matrix &operator=(const double** array);

        /**
         * @brief Adds two matrices.
         * @param other Matrix to add.
         * @return Resulting matrix after addition.
         */
        Matrix operator+(const Matrix &other) const;

        /**
         * @brief Subtracts another matrix from the current matrix.
         * @param other Matrix to subtract.
         * @return Resulting matrix after subtraction.
         */
        Matrix operator-(const Matrix &other) const;

        /**
         * @brief Multiplies two matrices.
         * @param other Matrix to multiply with.
         * @return Resulting matrix after multiplication.
         */
        Matrix operator*(const Matrix &other) const;

        /**
         * @brief Adds another matrix to the current matrix (in-place).
         * @param other Matrix to add.
         * @return Reference to the current matrix.
         */
        const Matrix &operator+=(const Matrix &other);

        /**
         * @brief Subtracts another matrix from the current matrix (in-place).
         * @param other Matrix to subtract.
         * @return Reference to the current matrix.
         */
        const Matrix &operator-=(const Matrix &other);

        /**
         * @brief Checks if two matrices are equal.
         * @param other Matrix to compare with.
         * @return True if matrices are equal, false otherwise.
         */
        bool operator==(const Matrix &other) const;

        /**
         * @brief Checks if two matrices are not equal.
         * @param other Matrix to compare with.
         * @return True if matrices are not equal, false otherwise.
         */
        bool operator!=(const Matrix &other) const;

        /**
         * @brief Checks if the matrix is a zero matrix.
         * @return True if the matrix is a zero matrix, false otherwise.
         */
        bool operator!() const;

        /**
         * @brief Computes the transpose of the matrix.
         * @return Transposed matrix.
         */
        Matrix transpose() const;

        /**
         * @brief Computes the row echelon form (REF) of the matrix.
         * @return Matrix in REF.
         */
        Matrix ref() const;

        /**
         * @brief Computes the reduced row echelon form (RREF) of the matrix.
         * @return Matrix in RREF.
         */
        Matrix rref() const;

        /**
        * @brief Creates a submatrix by removing the specified row and column from the current matrix.
        * @param row The index of the row to be removed (1-based index).
        * @param col The index of the column to be removed (1-based index).
        * @return Matrix A new matrix with the specified row and column removed.
        */
        Matrix operator()(int row, int col) const;

        /**
         * @brief Accesses an element of the matrix (modifiable).
         * @param row Row index.
         * @return Pointer to the element at the specified row.
         */
        double *operator[](int row);

        /**
         * @brief Accesses an element of the matrix (read-only).
         * @param row Row index.
         * @return Constant pointer to the element at the specified row.
         */
        const double *operator[](int row) const;

        /**
         * @brief Converts the matrix to a 2D array.
         * @return Pointer to the 2D array representation of the matrix.
         */
        operator double**() const;

        /**
         * @brief Multiplies a scalar with a matrix.
         * @param scalar Scalar value to multiply.
         * @param matrix Matrix to multiply with.
         * @return Resulting matrix after scalar multiplication.
         */
        friend Matrix operator*(double scalar, const Matrix &matrix);

        /**
         * @brief Outputs the matrix to an output stream.
         * @param os Output stream to write to.
         * @param matrix Matrix to output.
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Matrix &matrix);

        /**
         * @brief Inputs the matrix from an input stream.
         * @param is Input stream to read from.
         * @param matrix Matrix to populate.
         * @return Reference to the input stream.
         */
        friend std::istream& operator>>(std::istream& is, Matrix &matrix);
};

class SquareMatrix: private Matrix {
    public:
        SquareMatrix(int = 1);
        SquareMatrix(const SquareMatrix&);
        SquareMatrix(double**, int);
        ~SquareMatrix();

        double trace() const;
        double determinant() const;
        SquareMatrix inverse() const;
        SquareMatrix adjoint() const;
        bool isLTriangular() const;
        bool isUTriangular() const;
        bool isDiagonal() const;
        bool isSymmetric() const;

    friend Matrix operator*(double, const Matrix&);
    friend std::ostream& operator<<(std::ostream&, const Matrix&);
    friend std::istream& operator>>(std::istream&, Matrix&);
};

#endif