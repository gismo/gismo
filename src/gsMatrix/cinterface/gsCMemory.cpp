#include <gsCore/gsLinearAlgebra.h>
#include <gsMatrix/cinterface/gsCMemory.h>

// Copy From colwise pointer to gsMatrix
//mat = gismo::gsAsMatrix<>(data,rows,cols);
// Copy From rowwise pointer to gsMatrix
// mat = gismo::gsAsMatrix<>(data,rows,cols).transpose();

// copy std::vector
//double* data = new double[vector.size()];
//std::memcpy(data,vector.data(), vector.size());

void copy_matrix_c_array(const gismo::gsMatrix<>& matrix,
                 double* result,
                 int len_result)
{
    // copy matrix with check on amount of result elements provided

    for (int col = 0; col != matrix.cols(); col++)
    {
        for (int row = 0; row != matrix.rows(); row++)
        {
            int i = col * matrix.rows() + row;
            if (i < len_result)
            {
                result[i] = matrix(row, col);
            }
        }
    }
}

double* make_c_array(const gismo::gsMatrix<>& matrix)
{
    // code below is equivalent to
    // double* result = new double[matrix.size()];
    // copyRange(matrix.data(), result, matrix.size());
    //
    // I prefer to use the code below, because it is independent
    // of row or column mode of Eigen.

    double* result = new double[matrix.rows() * matrix.cols()];
    for (int col = 0; col != matrix.cols(); col++)
    {
        for (int row = 0; row != matrix.rows(); row++)
        {
            result[col * matrix.rows() + row] = matrix(row, col);
        }
    }

    return result;
}
