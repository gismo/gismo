
#ifndef GCMEMORY_H
#define GCMEMORY_H

void copy_matrix_c_array(const gismo::gsMatrix<>& matrix,
                         double* result,
                         int len_result);

double* make_c_array(const gismo::gsMatrix<double>& matrix);

#endif // GCMEMORY_H
