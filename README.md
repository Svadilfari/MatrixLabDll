# MatrixDll

This is a C++ library for **matrix operations** and realisation of **PCA algorithm**

## class Matrix:

### - get_rows 
Returns the number of rows in matrix
### - get_columns 
Returns the number of columns in matrix
### - get_matrix 
Returns the std::vector of std::vectors of matrix rows
### - get_row_i 
Get the number i
Returns Matrix -- vector of the elements of i_th row 
### - get_column_j
Get the number j
Returns Matrix -- vector of the elements of j_th column
### - is_square
Returns true if matrix is square
### - is_same_size
Get the Matrix B
Returns true if this matrix is the same size as B
### - can_be_mult 
Get the Matrix B
Returns true if it is possible to multiply this matrix and matrix B
### - is_vector
Returns true if this matrix is a vector
### - concatenate
Get the Matrix B
Returns the Matrix -- a concatenation of this matrix and B 
### - Det
Returns the determinant of the matrix
### - Inverse 
Returns the Matrix -- inversion of this matrix
### - Adamar
Get the Matrix B
Returns the Matrix -- Adamar multiplication of this matrix and B
### - Trace 
Returns the trace of this matrix
### - Norm
Returns the Norm of this matrix (if this matrix is a vector, returns the norm of vector)
### - MaxNorm
Returns the max element of this matrix
### - Rank 
Returns the rank of this matrix
### - Transpose
Returns the Matrix -- transposed version of this matrix

## Heirs of class Matrix:
### class IdMatriix
The identity matrix
### class UpperTriMatrix
The upper triangular matrix
### class LowerTriMatrix
The lower triangular matrix
### class DiagonalMatrix
The diagonal matrix
### class SymmetricMatrix
The symmetric matrix
