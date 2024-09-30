///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLBFGS                                                                    //
// http://www.loria.fr/~liuyang/software/HLBFGS/							 //
//                                                                           //
// HLBFGS is a hybrid L-BFGS optimization framework which unifies L-BFGS     //
// method, Preconditioned L-BFGS method and                                  //
// Preconditioned Conjugate Gradient method.                                 //
//                                                                           //
// Version 1.2                                                               //
// March 09, 2010                                                            //
//                                                                           //
// Copyright (C) 2009--2010                                                  //
// Yang Liu                                                                  //
//																			 //
// xueyuhanlang@gmail.com                                                    //
//                                                                           //
// HLBFGS is HLBFGS is freely available for non-commercial purposes.		 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Lite_Sparse_Matrix.h"

std::ostream & operator<<(std::ostream & s, Lite_Sparse_Matrix * A)
{
	s.precision(16);
	if (A == NULL)
	{
		s << "the matrix does not exist !\n ";
	}

	const int row = A->rows();
	const int col = A->cols();
	const int nonzero = A->get_nonzero();
	const int *rowind = A->get_rowind();
	const int *colptr = A->get_colptr();
	const double *values = A->get_values();

	s << "row :" << row << " col :" << col << " Nonzero: " << nonzero << "\n\n";

	s << "matrix --- (i, j, value)\n\n";

	SPARSE_STORAGE s_store = A->storage();
	int inc = (A->get_arraytype() == FORTRAN_TYPE ? -1 : 0);
	if (s_store == CCS)
	{
		int k = 0;
		for (int i = 1; i < col + 1; i++)
		{
			for (int j = 0; j < colptr[i] - colptr[i - 1]; j++)
			{
				s << rowind[k] + inc << " " << i - 1 << " " << std::scientific
					<< values[k] << "\n";
				k++;
			}
		}
	}
	else if (s_store == CRS)
	{
		int k = 0;
		for (int i = 1; i < row + 1; i++)
		{
			for (int j = 0; j < rowind[i] - rowind[i - 1]; j++)
			{
				s << i - 1 << " " << colptr[k] + inc << " " << std::scientific
					<< values[k] << "\n";
				k++;
			}
		}
	}
	else if (s_store == TRIPLE)
	{
		for (int k = 0; k < nonzero; k++)
		{
			s << rowind[k] + inc << " " << colptr[k] + inc << " "
				<< std::scientific << values[k] << "\n";
		}
	}

	return s;
}
