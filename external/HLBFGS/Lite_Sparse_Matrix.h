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


#ifndef LITE_SPARSE_MATRIX_H
#define LITE_SPARSE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cassert>
#include <ostream>
#include "Sparse_Entry.h"

// \addtogroup MathSuite
//@{

//! Symmetric status
enum SYMMETRIC_STATE
{
	NOSYM, /*!< general case */
	SYM_UPPER, /*!< symmetric (store upper triangular part) */
	SYM_LOWER, /*!< symmetric (store lower triangular part) */
	SYM_BOTH
	/*!< symmetric (store both upper and lower triangular part) */
};

//!     Storage
enum SPARSE_STORAGE
{
	CCS, /*!< compress column format */
	CRS, /*!< compress row format */
	TRIPLE
	/*!< row-wise coordinate format */
};

//! Array type
enum ARRAYTYPE
{
	FORTRAN_TYPE, /*!< the index starts from 1 */
	C_TYPE
	/*!< the index starts from 0 */
};

//////////////////////////////////////////////////////////////////////////
//! Lite Sparse Matrix Class
class Lite_Sparse_Matrix
{
private:

	//!     Status for creating sparse solver
	enum STORE_STATE
	{
		ENABLE, DISABLE, LOCK
	};

	STORE_STATE state_fill_entry;
	SYMMETRIC_STATE sym_state;
	SPARSE_STORAGE s_store;
	ARRAYTYPE arraytype;

	int nrows; //!< number of rows
	int ncols; //!< number of columns
	int nonzero; //!< number of nonzeros
	//! pointers to where columns begin in rowind and values 0-based, length is (col+1)
	/*!
	* When s_store is CRS, colptr stores column indices;
	*/
	std::vector<int> colptr;
	//! row indices, 0-based
	/*!
	* When s_store is CRS, rowind stores row-pointers
	*/
	std::vector<int> rowind;
	std::vector<double> values; //!< nonzero values of the sparse matrix

	std::vector<std::vector<Sparse_Entry> > entryset; //!< store temporary sparse entries

	std::vector<double> diag; //! special usage for some libraries

	bool save_diag_separetely;

public:

	//! Sparse matrix constructor
	/*!
	* \param m row dimension
	* \param n column dimension
	* \param symmetric_state
	* \param m_store the storage format
	* \param atype Fortran or C type of array
	*/
	Lite_Sparse_Matrix(int m, int n, SYMMETRIC_STATE symmetric_state = NOSYM,
		SPARSE_STORAGE m_store = CCS, ARRAYTYPE atype = C_TYPE,
		bool save_diag = false) :
	nrows(m), ncols(n), sym_state(symmetric_state), s_store(m_store),
		arraytype(atype), save_diag_separetely(save_diag), nonzero(0),
		state_fill_entry(DISABLE)
	{
		if (m != n)
		{
			symmetric_state = NOSYM;
		}

		int nn = (m_store == CCS ? ncols : nrows);
		entryset.resize(nn);
		if (save_diag_separetely)
		{
			diag.resize(nrows < ncols ? nrows : ncols);
			std::fill(diag.begin(), diag.end(), 0.0);
		}
	}

	//! Sparse matrix destructor
	~Lite_Sparse_Matrix()
	{
		clear_mem();
	}
	//! Start to build sparse matrix pattern
	inline void begin_fill_entry()
	{
		state_fill_entry = ENABLE;
	}
	//! Construct sparse pattern
	void end_fill_entry()
	{
		assert (state_fill_entry == ENABLE);

		clear_mem();

		state_fill_entry = LOCK;

		int inc = (arraytype == FORTRAN_TYPE ? 1 : 0);

		if (s_store == CCS)
		{
			//construct map and ccs matrix
			int i, j, k = 0;
			colptr.resize(ncols + 1);
			colptr[0] = inc;
			for (j = 1; j < ncols + 1; j++)
			{
				colptr[j] = (int) entryset[j - 1].size() + colptr[j - 1];
			}

			nonzero = colptr[ncols];

			if (nonzero > 0)
			{
				rowind.resize(nonzero);
				values.resize(nonzero);

				for (j = 0; j < ncols; j++)
				{
					for (i = 0; i < colptr[j + 1] - colptr[j]; i++)
					{
						rowind[k] = entryset[j][i].index + inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}
			}
		}
		else if (s_store == CRS)
		{
			//construct map and crs matrix
			int i, j, k = 0;
			rowind.resize(nrows + 1);
			rowind[0] = inc;
			for (j = 1; j < nrows + 1; j++)
			{
				rowind[j] = (int) entryset[j - 1].size() + rowind[j - 1];
			}
			nonzero = rowind[nrows];
			if (nonzero > 0)
			{
				colptr.resize(nonzero);
				values.resize(nonzero);

				for (j = 0; j < nrows; j++)
				{
					for (i = 0; i < rowind[j + 1] - rowind[j]; i++)
					{
						colptr[k] = entryset[j][i].index + inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}
			}
		}
		else if (s_store == TRIPLE)
		{
			int i, j, k = 0;
			nonzero = 0;
			for (i = 0; i < nrows; i++)
			{
				nonzero += (int) entryset[i].size();
			}

			if (nonzero > 0)
			{
				rowind.resize(nonzero);
				colptr.resize(nonzero);
				values.resize(nonzero);

				for (i = 0; i < nrows; i++)
				{
					int jsize = (int) entryset[i].size();
					for (j = 0; j < jsize; j++)
					{
						rowind[k] = i + inc;
						colptr[k] = entryset[i][j].index + inc;
						values[k] = entryset[i][j].value;
						k++;
					}
				}
			}
		}
		entryset.clear();
	}

	//! Fill matrix entry \f$  Mat_{row_index, col_index} += val \f$
	void fill_entry(int row_index, int col_index, double val = 0)
	{
		if (row_index >= nrows || col_index >= ncols)
			return;

		if (save_diag_separetely && row_index == col_index)
		{
			diag[row_index] += val;
		}

		if (sym_state == NOSYM)
		{
			fill_entry_internal(row_index, col_index, val);
		}
		else if (sym_state == SYM_UPPER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(row_index, col_index, val);
			}
			else
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
		else if (sym_state == SYM_LOWER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
			else
			{
				fill_entry_internal(row_index, col_index, val);
			}
		}
		else if (sym_state == SYM_BOTH)
		{
			fill_entry_internal(row_index, col_index, val);

			if (row_index != col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
	}

	//fill the diagonal entry
	inline void fill_diag(int diagid, double v = 0)
	{

		if (diag.size() == 0)
		{
			diag.resize(nrows < ncols ? nrows : ncols);
			std::fill(diag.begin(), diag.end(), 0.0);
		}
		diag[diagid] += v;
	}

	//! get the number of nonzeros
	inline int get_nonzero()
	{
		return nonzero;
	}
	//! get the row dimension
	inline int rows()
	{
		return nrows;
	}
	//! get the column dimension
	inline int cols()
	{
		return ncols;
	}
	//! return the symmetric state
	inline bool issymmetric()
	{
		return sym_state != NOSYM;
	}

	//! tell whether the matrix is upper or lower symmetric
	inline bool issym_store_upper_or_lower()
	{
		return (sym_state == SYM_LOWER) || (sym_state == SYM_UPPER);
	}

	//! return symmetric state
	inline SYMMETRIC_STATE symmetric_state()
	{
		return sym_state;
	}

	//! tell whether the matrix is square
	inline bool issquare()
	{
		return nrows == ncols;
	}

	//! return the storage format
	inline SPARSE_STORAGE storage()
	{
		return s_store;
	}

	//! return array type
	inline ARRAYTYPE get_arraytype()
	{
		return arraytype;
	}

	//! get rowind
	inline int *get_rowind()
	{
		return &rowind[0];
	}

	inline const int *get_rowind() const
	{
		return &rowind[0];
	}

	//! get colptr
	inline int *get_colptr()
	{
		return &colptr[0];
	}

	inline const int *get_colptr() const
	{
		return &colptr[0];
	}

	//! get the values array
	inline double *get_values()
	{
		return &values[0];
	}

	inline const double *get_values() const
	{
		return &values[0];
	}

	//! get the diagonal array
	inline double *get_diag()
	{
		return &diag[0];
	}

	inline const double *get_diag() const
	{
		return &diag[0];
	}

	//////////////////////////////////////////////////////////////////////////
private:
	//! Clear memory
	inline void clear_mem()
	{
		colptr.clear();
		rowind.clear();
		values.clear();
	}

	//! fill matrix entry (internal) \f$ Mat[rowid][colid] += val \f$
	bool fill_entry_internal(int row_index, int col_index, double val = 0)
	{
		assert (state_fill_entry == ENABLE);

		int search_index = (s_store == CCS ? row_index : col_index);
		int pos_index = (s_store == CCS ? col_index : row_index);

		Sparse_Entry forcompare(search_index);

		std::vector<Sparse_Entry>::iterator iter = std::lower_bound(
			entryset[pos_index].begin(), entryset[pos_index].end(),
			forcompare);
		if (iter != entryset[pos_index].end())
		{
			if (iter->index == search_index)
			{
				iter->value += val;
			}
			else
				entryset[pos_index].insert(iter,
				Sparse_Entry(search_index, val));
		}
		else
		{
			entryset[pos_index].push_back(Sparse_Entry(search_index, val));
		}
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
};

//! print sparse matrix
std::ostream & operator<<(std::ostream & s, const Lite_Sparse_Matrix * A);

//@}


#endif //Lite_Sparse_Matrix_H
