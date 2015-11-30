/** @file gsAsMatrix.h

    @brief Wraps pointers as matrix objects

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

# pragma once

// Assumes that Eigen library has been already included

namespace gismo
{

template<class T, int _Rows, int _Cols>
class gsAsConstMatrix;

/** \brief Creates a mapped object or data pointer to a matrix without
    copying data.
    
   This allows for re-indexing the matrix. No copying is taking place
   and the original matrix remains untached.

   \tparam T coefficient type
   \ingroup Matrix
*/
template<class T, int _Rows, int _Cols>
class gsAsMatrix : public Eigen::Map< Eigen::Matrix<T,_Rows,_Cols> >
{
public:
    typedef Eigen::Map< Eigen::Matrix<T,_Rows,_Cols> > Base;

    // type of row minor matrix: rows reduced by one
    typedef gsMatrix< T, ChangeDim<_Rows, -1>::D, _Cols>
        RowMinorMatrixType;

    // type of col minor matrix: cols reduced by one
    typedef gsMatrix< T, _Rows, ChangeDim<_Cols, -1>::D>
        ColMinorMatrixType;

public:
    gsAsMatrix( std::vector<T> & v, index_t n, index_t m)
    : Base( &v[0], n, m)
    { 
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
        GISMO_ASSERT( m*n <= index_t(v.size()), "Not enough coefficients in vector to map." ); 
    }

    gsAsMatrix( std::vector<T> & v)
    : Base( &v[0], 1, v.size() ) 
    {  
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
    }

    gsAsMatrix( T * pt, unsigned n, unsigned m)
    : Base( pt, n, m) {  }

#ifdef _MSC_VER
    template <class EigenExpr>
    gsAsMatrix& operator= (const EigenExpr & other) 
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif

    /// Returns the ith row minor, i.e. the matrix after removing row
    /// \a i from the matrix. After the operation the row size of the
    /// matrix is one less.
    void rowMinor(index_t i, RowMinorMatrixType & result ) const
    {
        const index_t mrows = this->rows()-1;
        GISMO_ASSERT( i <= mrows, "Invalid row." );
        result.resize(mrows, Eigen::NoChange);
        result.topRows(i)          = this->topRows(i);
        result.bottomRows(mrows-i) = this->bottomRows(mrows-i);
    }
    
    /// Returns the jth column minor, i.e. the matrix after removing row
    /// \a j from the matrix. After the operation the column size of the
    /// matrix is one less.
    void colMinor(index_t j, ColMinorMatrixType & result ) const
    {
        const index_t mcols = this->cols()-1;
        GISMO_ASSERT( j <= mcols, "Invalid column." );
        result.resize(Eigen::NoChange, mcols);
        result.leftCols(j)        = this->leftCols(j);
        result.rightCols(mcols-j) = this->rightCols(mcols-j);
    }

    operator gsAsConstMatrix<T,_Rows,_Cols>() const {return gsAsConstMatrix<T,_Rows,_Cols>(&(*this)(0,0),Base::rows(),Base::cols());}
private:
    gsAsMatrix() { }
};

/** \brief Creates a mapped object or data pointer to a const matrix without
    copying data.
    
   This allows for re-indexing the matrix. No copying is taking place
   and the original matrix remains untached.

   \tparam T coefficient type
   \ingroup Matrix
*/
template<class T, int _Rows, int _Cols>
class gsAsConstMatrix : public Eigen::Map< const Eigen::Matrix<T,_Rows,_Cols> >
{
public:
    typedef Eigen::Map<const Eigen::Matrix<T,_Rows,_Cols> > Base;

public:

    gsAsConstMatrix( const std::vector<T> & v, index_t n, index_t m)
    : Base( &v[0], n, m)
    { 
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
        GISMO_ASSERT( m*n <= index_t(v.size()), "Not enough coefficients in vector to map." ); 
    }

    gsAsConstMatrix( const std::vector<T> & v)
    : Base( &v[0], 1, v.size() ) 
    {  
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
    }

    gsAsConstMatrix( const T * pt, unsigned n, unsigned m)
    : Base( pt, n, m) {  }

private:
    gsAsConstMatrix() { }
};

/** \brief Creates a mapped object or data pointer to a vector without
    copying data.
    
   This allows for re-indexing the matrix. No copying is taking place
   and the original matrix remains untached.

   \tparam T coefficient type
   \ingroup Matrix
 */
template<class T, int _Rows=Dynamic>
class gsAsVector : public Eigen::Map< Eigen::Matrix<T,_Rows,1> >
{
public:
    typedef Eigen::Map< Eigen::Matrix<T,_Rows,1> > Base;

    // Type for treating a vector as a permutation matrix
    typedef Eigen::PermutationMatrix<_Rows> Permutation;

public:
    gsAsVector( std::vector<T> & v)
    : Base( &v[0], v.size(), 1 ) 
    {  
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
    }

    gsAsVector( T * pt, unsigned n)
    : Base( pt, n) {  }

#ifdef _MSC_VER
    template <class EigenExpr>
    gsAsVector& operator= (const EigenExpr & other) 
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif

private:
    gsAsVector() { }
};

/** \brief Creates a mapped object or data pointer to a const vector without
    copying data.
    
   This allows for re-indexing the matrix. No copying is taking place
   and the original matrix remains untached.

   \tparam T coefficient type
   \ingroup Matrix
 */
template<class T, int _Rows=Dynamic>
class gsAsConstVector : public Eigen::Map< const Eigen::Matrix<T,_Rows,1> >
{
public:
    typedef Eigen::Map<const Eigen::Matrix<T,_Rows,1> > Base;

public:

    gsAsConstVector( const std::vector<T> & v)
    : Base( &v[0], v.size(), 1) 
    {  
        GISMO_ASSERT( v.size() != 0, "Tried to map an empty vector." ); 
    }

    gsAsConstVector( const T * pt, unsigned n)
    : Base( pt, n, 1) {  }

private:
    gsAsConstVector() { }
};


/// Utility to make a matrix out of an iterator to values
template<class T, class iterator>
typename gsMatrix<T>::uPtr makeMatrix(iterator it, index_t n, index_t m)
{
    typename gsMatrix<T>::uPtr result ( new gsMatrix<T>(n,m) );
    for ( index_t i = 0; i!=n; ++i)
        for ( index_t j = 0; j!=m; ++j)
            (*result)(i,j)= *(it++);
    return result;
}


} // namespace gismo
