/** @file gsMatrix.h

    @brief Provides declaration of Matrix class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

# pragma once

namespace gismo
{

/** @brief
    A matrix with arbitrary coefficient type and fixed or dynamic size.

    This class provides an interface to Eigen::Matrix from the Eigen
    linear algebra library. Most operations from Eigen are supported
    on a gsMatrix. 

    See therefore also the Eigen documentation for dense matrices,
    http://eigen.tuxfamily.org/dox/group__QuickRefPage.html

    \tparam T coefficient type
    \tparam _Rows number of rows: an integer or \c Dynamic
    \tparam _Cols number of rows: an integer or \c Dynamic
    \tparam _Options further options; see Eigen documentation

    \ingroup Matrix
*/
template<class T, int _Rows, int _Cols, int _Options>
class gsMatrix : public Eigen::Matrix<T,_Rows, _Cols, _Options>
{
public:
    // Base is the dense matrix class of Eigen
    typedef Eigen::Matrix<T,_Rows, _Cols, _Options> Base;

    // The type of the coefficients of the matrix
    typedef T Scalar_t;

    // Type pointing to a block view of the matrix
    typedef gsMatrixBlockView<Base> BlockView;

    // Type pointing to a block of the matrix
    typedef Eigen::Block<Base> Block;

    // Type pointing to a (const) block of the matrix
    typedef Eigen::Block<const Base> constBlock;
    
    // Type pointing to a row of the matrix
    typedef Eigen::Block<Base, 1, _Cols, false> Row;

    // Type pointing to a (const) row of the matrix
    typedef Eigen::Block<const Base, 1, _Cols, false> constRow;

    // Type pointing to a set of successive rows of the matrix
    typedef Eigen::Block<Base, Dynamic, _Cols, false> Rows;

    // Type pointing to a a set of successive (const) rows of the matrix
    typedef Eigen::Block<const Base, Dynamic, _Cols, false> constRows;

    // Type pointing to a column of the matrix
    typedef Eigen::Block<Base, _Rows, 1, true > Column;

    // Type pointing to a (const) column of the matrix
    typedef Eigen::Block<const Base, _Rows, 1, true > constColumn;

    // Type pointing to a set of successive columns of the matrix
    typedef Eigen::Block<Base, _Rows, Dynamic, true > Columns;

    // Type pointing to a set of successive (const) columns of the matrix
    typedef Eigen::Block<const Base, _Rows, Dynamic, true > constColumns;

    // Type pointing to the transpose of the matrix
    typedef Eigen::Transpose<Base> Tr;

    // Type pointing to the (const) transpose of the matrix
    typedef const Eigen::Transpose<const Base> constTr;

    // Type refering to any possible Eigen type that can be copied
    // into a gsMatrix
    typedef Eigen::Ref<Base> Ref;
    
    // Type refering to any (const) possible Eigen types that can be
    // copied into a gsMatrix
    typedef const Eigen::Ref<const Base> constRef;

    /// Shared pointer for gsMatrix
    typedef memory::shared_ptr< gsMatrix > Ptr;

    /// Unique pointer for gsMatrix
//  #ifdef __GXX_EXPERIMENTAL_CXX0X__
//     typedef std::unique_ptr< gsMatrix > uPtr;
//  #else
    typedef memory::auto_ptr< gsMatrix > uPtr;
// #endif

    // type of first minor matrix: rows and cols reduced by one
    typedef gsMatrix< T, ChangeDim<_Rows, -1>::D, ChangeDim<_Cols, -1>::D, _Options>
        FirstMinorMatrixType;

    // type of row minor matrix: rows reduced by one
    typedef gsMatrix< T, ChangeDim<_Rows, -1>::D, _Cols, _Options>
        RowMinorMatrixType;

    // type of col minor matrix: cols reduced by one
    typedef gsMatrix< T, _Rows, ChangeDim<_Cols, -1>::D, _Options>
        ColMinorMatrixType;

public:

    gsMatrix() ;
    gsMatrix(const Base& a) ;
    gsMatrix(int rows, int cols) ;

    /// This constructor allows constructing a gsMatrix from Eigen expressions
    template<typename OtherDerived>
    gsMatrix(const Eigen::EigenBase<OtherDerived>& other) : Base(other) { }

    /// This constructor allows constructing a gsMatrix from Eigen expressions
    template<typename OtherDerived>
    gsMatrix(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) { }

    /// This constructor allows constructing a gsMatrix from Eigen expressions
    template<typename OtherDerived>
    gsMatrix(const Eigen::ReturnByValue<OtherDerived>& other) : Base(other) { }

    /// move constructor
    gsMatrix(gsMovable< gsMatrix > other)
    {
        this->swap( other.ref() );
    }

    /// constructor by swapping a unique pointer
    gsMatrix(uPtr const & other)
    {
        // NOTE: due to a language edge case, we can not pass the uPtr by value
        // into this constructor as we should. Assuming that uPtr is an auto_ptr,
        // this fails at least on gcc due to implicit conversion rules (no
        // two user-defined conversions in an initialization).
        this->swap( *other );
    }

    ~gsMatrix() ;

    // Using the assignment operators of Eigen
    // Note: using Base::operator=; is ambiguous in MSVC
#ifdef _MSC_VER // && !__INTEL_COMPILER
    template <class EigenExpr>
    gsMatrix& operator= (const EigenExpr & other) 
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif

    /// This method allows to swap with another matrix
    gsMatrix& operator= (gsMovable< gsMatrix > other)
    {
        this->resize(0,0);
        this->swap( other.ref() );
        return *this;
    }

    gsMatrix& operator= (uPtr other)
    {
        this->resize(0,0);
        this->swap( *other );
        return *this;
    }

    /// \brief Returns the \a i-th element of the vectorization of the matrix
    T   at (index_t i) const { return *(this->data()+i);}

    /// \brief Returns the \a i-th element of the vectorization of the matrix
    T & at (index_t i)       { return *(this->data()+i);}

    // \brief Returns the last element of the matrix (maximum row and column)
    //T   lastCoeff() { return *(this->data()+this->size()-1);}

    /// \brief Returns the matrix resized to n x m matrix (data is not copied)
    gsAsMatrix<T, Dynamic, Dynamic> reshape(index_t n, index_t m )
    { return gsAsMatrix<T, Dynamic, Dynamic>(this->data(), n, m); }

    /// \brief Returns column \a c of the matrix resized to n x m matrix
    gsAsMatrix<T, Dynamic, Dynamic> reshapeCol( index_t c, index_t n, index_t m )
    { return gsAsMatrix<T, Dynamic, Dynamic>(this->col(c).data(), n, m); }

    /// \brief Returns column \a c of the matrix resized to n x m matrix
    gsAsConstMatrix<T, Dynamic, Dynamic> reshapeCol( index_t c, index_t n, index_t m ) const
    { return gsAsConstMatrix<T, Dynamic, Dynamic>(this->col(c).data(), n, m); }

    /// \brief Returns the entries of the matrix resized to a n*m vector column-wise
    gsAsVector<T, Dynamic> asVector()
    { return gsAsVector<T, Dynamic>(this->data(), this->rows()*this->cols() ); }

    /// \brief Returns the entries of the matrix resized to a (const) n*m vector column-wise
    gsAsConstVector<T, Dynamic> asVector() const
    { return gsAsConstVector<T, Dynamic>(this->data(), this->rows()*this->cols() ); }

    /// Returns the entries of the matrix resized to a 1 x n*m
    /// row-vector column-wise
    gsAsMatrix<T, 1, Dynamic> asRowVector()
    { return gsAsMatrix<T, 1, Dynamic>(this->data(), 1, this->rows()*this->cols() ); }

    /// Returns the entries of the matrix resized to a (const) 1 x n*m
    /// row-vector column-wise
    gsAsConstMatrix<T, 1, Dynamic> asRowVector() const
    { return gsAsConstMatrix<T, 1, Dynamic>(this->data(), 1, this->rows()*this->cols() ); }

    /// Get a submatrix consisting of the columns indexed by the vector cols
    gsMatrix<T,_Rows,Dynamic> * submatrixCol( std::vector<index_t> & cols ) const
    {
        const index_t cc= cols.size();
        gsMatrix<T,_Rows,Dynamic> res = new gsMatrix<T,_Rows,Dynamic>(_Rows, cc);
        for ( index_t i = 0; i!= cc; ++i )
            res->col(i) = this->col(i);

        return res;
    }

    /// Removes column \a i from the matrix. After the operation the
    /// column size of the matrix is one less.
    void removeCol( index_t i )
    {
        index_t cc= this->cols();
        GISMO_ASSERT( i < cc, "Invalid column." );
        for ( index_t c = i+1; c!= cc; ++c )
            this->col(c-1) = this->col(c);
        this->conservativeResize(this->rows(), cc-1);
    }

    /// Returns the (i,j)-minor, i.e. the matrix after removing row
    /// \a i and column \a j from the matrix. After the operation the
    /// row and column size of the matrix is one less.
    void firstMinor(index_t i, index_t j, FirstMinorMatrixType & result ) const
    {
        const index_t mrows = this->rows()-1, 
            mcols = this->cols()-1;
        GISMO_ASSERT( i <= mrows, "Invalid row." );
        GISMO_ASSERT( j <= mcols, "Invalid column." );
        result.resize(mrows,mcols);
        result.block(0,0,i,j)             = this->block(0,0,i,j);
        result.block(i,0,mrows-i,j)       = this->block(i+1,0,mrows-i,j);
        result.block(0,j,i,mcols-j)       = this->block(0,j+1,i,mcols-j);
        result.block(i,j,mrows-i,mcols-j) = this->block(i+1,j+1,mrows-i,mcols-j);
    }

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

    void duplicateRow( index_t k )
    {
        this->conservativeResize(this->rows() + 1, this->cols()); 

        /*
        // Test this
        this->bottomRows(this->rows() - k ) = 
        this->middleRows(this->rows() - k, k+1 );
        
        this->row(k+1) = this->row(k);
        return;

        //*/

        for (index_t i = this->rows() - 1; i > k+1 ; --i)
            this->row(i).swap(this->row(i-1));

        this->row(k+1) = this->row(k);
    }

    /// Clone function. Used to make a copy of the matrix
    gsMatrix * clone() const;

    /// Return a block view of the matrix with \a rowSizes and \a colSizes
    BlockView blockView(const gsVector<index_t> & rowSizes, 
                        const gsVector<index_t> & colSizes)
    {
        return BlockView(*this, rowSizes, colSizes);
    }

    /// Sorts rows of matrix by column \em j.
    void sortByColumn( index_t j )
    {
        GISMO_ASSERT( j < this->cols(), "Invalid column.");

        unsigned lastSwapDone = this->rows() - 1;
        unsigned lastCheckIdx = lastSwapDone;

        bool didSwap;
        gsMatrix<T> tmp(1, this->cols() );
        do{
            didSwap = false;
            lastCheckIdx = lastSwapDone;

            for( unsigned i=0; i < lastCheckIdx; i++)
                if( this->coeff(i,j) > this->coeff(i+1,j) )
                {
                    tmp.row(0) = this->row(i);
                    this->row(i) = this->row(i+1);
                    this->row(i+1) = tmp.row(0);

                    didSwap = true;
                    lastSwapDone = i;
                }
        }while( didSwap );
    }

/// \brief Transposes in place the matrix block-wise. The matrix is
//  treated a rows() x (cols()/colBlock) block matrix, and every block
//  of size rows() x colBlock is transposed in place
void blockTransposeInPlace(const index_t colBlock)
{
    const index_t nc = this->cols();
    const index_t nr = this->rows();
    
    GISMO_ASSERT( nc % colBlock == 0,
                  "The blocksize is not compatible with number of columns.");

    if (nr == 1 || colBlock == 1)
    {
        this->resize(colBlock, this->size()/colBlock);
    }
    else if ( nr == colBlock )
    {
        for (index_t j = 0; j!= nc; j+=colBlock)
            this->middleCols(j,colBlock).template triangularView<Eigen::StrictlyUpper>()
                .swap( this->middleCols(j,colBlock).transpose() );
    }
    else
    {
        Eigen::Map<Base> m(this->data(), nr, nc);
        this->resize(colBlock, this->size()/colBlock);
        
        index_t i = 0;
        for (index_t j = 0; j!= nc; j+=colBlock, i+=nr)
            this->middleCols(i,nr) = m.middleCols(j,colBlock).transpose().eval();
    }
}

}; // class gsMatrix



////////////////////////////////////////////////
////////////////////////////////////////////////




template<class T, int _Rows, int _Cols, int _Options> inline
gsMatrix<T,_Rows, _Cols, _Options>::gsMatrix() { }
 
template<class T, int _Rows, int _Cols, int _Options> inline
gsMatrix<T,_Rows, _Cols, _Options>::gsMatrix(const Base& a) : Base(a) { }

template<class T, int _Rows, int _Cols, int _Options> inline
gsMatrix<T,_Rows, _Cols, _Options>::gsMatrix(int rows, int cols) : Base(rows,cols) { }
    
// template<class T, int _Rows, int _Cols, int _Options>
//  template<typename OtherDerived> 
// gsMatrix<T,_Rows, _Cols, _Options>::gsMatrix(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) { }
    

template<class T, int _Rows, int _Cols, int _Options> inline
gsMatrix<T,_Rows, _Cols, _Options>::~gsMatrix() { }
  
/// Clone function. Used to make a copy of the matrix
template<class T, int _Rows, int _Cols, int _Options> inline
gsMatrix<T,_Rows, _Cols, _Options> * gsMatrix<T,_Rows, _Cols, _Options>::clone() const
{ return new gsMatrix<T,_Rows, _Cols, _Options>(*this); }
  


}; // namespace gismo
