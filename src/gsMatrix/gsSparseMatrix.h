/** @file gsSparseMatrix.h

    @brief Provides declaration of the gsSparseMatrix class.

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

/** @brief Collection of available sparse linear solvers

    Example of usage:
    \code
    gsSparseMatrix<> M; // sparse system matrix
    gsMatrix<> b; // right-hand side
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute(M);
    gsMatrix<> x = solver.solve(b);
    \endcode

    The template arguments are the same as the ones for gismo::gsSparseMatrix
   See also http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
*/
template<typename T=real_t, int _Options=0, typename _Index = index_t>
struct gsSparseSolver
{
    // Note: IncompleteILU is not compatible with
    // Eigen::ConjugateGradient because this preconditionner does not
    // preserve symmetry.

    /// Congugate gradient without preconditioner (identity as preconditioner) 
    /// The matrix is assumed symmetric and only the lower trinagular part is used
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<T,_Options,_Index>,
            Eigen::Lower, Eigen::IdentityPreconditioner> CGIdentity;

    /// Congugate gradient with diagonal (Jacobi) preconditioner
    /// The matrix is assumed symmetric and only the lower trinagular part is used
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<T,_Options,_Index>, 
            Eigen::Lower, Eigen::DiagonalPreconditioner<T> > CGDiagonal;

    /// BiCGSTAB with Incomplete LU factorization with dual-threshold strategy
    typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<T,_Options,_Index>,
                            Eigen::IncompleteLUT<T> > BiCGSTABILUT;

    /// Direct LDLt factorization
    typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<T,_Options,_Index> > SimplicialLDLT;

    /// Sparse LU solver
    typedef Eigen::SparseLU<Eigen::SparseMatrix<T,_Options,_Index>,
                            Eigen::COLAMDOrdering<_Index> > LU;

    /// Sparse QR solver
    typedef Eigen::SparseQR<Eigen::SparseMatrix<T,_Options,_Index>,
                            Eigen::COLAMDOrdering<_Index> > QR;
};

/**
   @brief Class that provides a container for triplets (i,j,value) to be
   filled in a sparse matrix.  

   Constructing a sparse matrix from triplets is much faster than
   inserting directly.  Use gsSparseMatrix().setFrom(gsSparseEntries)
   to pass the triplets to the matrix.

   \tparam T coefficient type
   \ingroup Matrix
*/
template<class T>
class gsSparseEntries : public std::vector<Eigen::Triplet<T,index_t> >
{
public:
    typedef Eigen::Triplet<T,index_t> Triplet;
    typedef std::vector<Eigen::Triplet<T,index_t> > Base;

public:
    gsSparseEntries() ;

    ~gsSparseEntries() ;

    inline void add( int i, int j, T value )
    { this->push_back( Triplet(i,j,value) ); }

};

/** @brief Sparse matrix class, based on Eigen::SparseMatrix.
 *
 * See http://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html
 * for Eigen's sparse matrix manipulations and
 * http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html for
 * documentation of the Eigen::SparseMatrix class.
 *
 * Remarks:
 *
 * An entry of the gsSparseMatrix can be accessed by <b>coeff( index
 * row, index col )</b> or just with the operator <b>( index row, index col )</b>.\n
 * An entry can be changed with either
 * <b>coeffRef( index row, index col)</b> or operator <b>( index row, index col )</b>.\n
 *
 *   \tparam T coefficient type
 *   \tparam _Option zero is ColMajor order.
 *   \tparam _Index index type
 *  \ingroup Matrix
 */

// Export the result to a file: saveAsBitmap(...);

template<typename T, int _Options, typename _Index>
class gsSparseMatrix : public Eigen::SparseMatrix<T,_Options,_Index>
{
public:
    typedef Eigen::SparseMatrix<T,_Options,_Index> Base;

    // Type pointing to a block of the sparse matrix
    typedef typename Eigen::Block<Base> Block;

    // Type pointing to a const block of the sparse matrix
    typedef typename Eigen::Block<const Base> constBlock;

    // Type pointing to a block view of the sparse matrix
    typedef gsMatrixBlockView<Base> BlockView;

    /// Shared pointer for gsSparseMatrix
    typedef memory::shared_ptr< gsSparseMatrix > Ptr;

    /// Type of the full view of the matrix, for the case when only
    /// the lower diagonal part is stored
    typedef typename Eigen::SparseSelfAdjointView<Base, Lower> fullView;

    /// Type of the full view of the matrix, for the case when only
    /// the lower diagonal part is stored
    typedef typename Eigen::SparseSelfAdjointView<const Base, Lower> constFullView;

public:
    gsSparseMatrix() ;

    gsSparseMatrix(_Index rows, _Index cols) ;

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const Eigen::EigenBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from a selfadjoint view
    template<typename OtherDerived, unsigned int UpLo>
    gsSparseMatrix(const Eigen::SparseSelfAdjointView<OtherDerived, UpLo>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const Eigen::MatrixBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseMatrix from another sparse expression
    template<typename OtherDerived> 
    gsSparseMatrix(const Eigen::SparseMatrixBase<OtherDerived>& other)  : Base(other) { } 

    /// This constructor allows constructing a gsSparseMatrix from Eigen expressions
    template<typename OtherDerived>
    gsSparseMatrix(const Eigen::ReturnByValue<OtherDerived>& other)  : Base(other) { }

    gsSparseMatrix(gsMovable< gsSparseMatrix > other)
    {
        this->swap( other.ref() );
    }

    // Using the assignment operators of Eigen
    // Note: using Base::operator=; is ambiguous in MSVC
#ifdef _MSC_VER
    template <class EigenExpr>
    gsSparseMatrix& operator= (const EigenExpr & other) 
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif

    /// This method allows to swap with another vector
    gsSparseMatrix& operator=(gsMovable< gsSparseMatrix > other)
    {
        this->resize(0,0);
        this->swap( other.ref() );
        return *this;
    }


    ~gsSparseMatrix() ;

    void setFrom( gsSparseEntries<T> const & entries) ;

    inline T   at (_Index i, _Index j = 0) const { return this->coeff(i,j); }
    inline T & at (_Index i, _Index j = 0) { return this->coeffRef(i,j); }

    inline T    operator () (_Index i, _Index j = 0) const { return this->coeff(i,j); }
    inline T  & operator () (_Index i, _Index j = 0) { return this->coeffRef(i,j); }

    /// Clone function. Used to make a copy of the matrix
    gsSparseMatrix * clone() const ;

    /// Return a block view of the matrix with \a rowSizes and \a colSizes
    BlockView blockView(const gsVector<index_t> & rowSizes, 
                        const gsVector<index_t> & colSizes)
    {
        return BlockView(*this, rowSizes, colSizes);
    }


    /// Returns a pointer wrapped as a gsAsConstVector, which contains
    /// the number of non-zero entries per column. Note that the
    /// matrix must be uncompressed format for this to work
    gsAsConstVector<_Index> nonZerosPerCol()
    {
        if ( this->isCompressed() )
        {
            gsWarn<<"nonZerosPerCol(): Uncompressing the gsSparseMatrix.\n";
            this->uncompress();
        }
        return gsAsConstVector<_Index>(this->innerNonZeroPtr(),this->innerSize());
    }

}; // class gsSparseMatrix



////////////////////////////////////////////////
////////////////////////////////////////////////




template<class T> inline
gsSparseEntries<T>::gsSparseEntries() : Base() { }

template<class T> inline
gsSparseEntries<T>::~gsSparseEntries() {  }

//template<class T>
//inline void gsSparseEntries<T>::add( int const& i, int const& j, const T & value)
//        { this->push_back( Triplet(i,j,value) ); }

//
// gsSparseMatrix //////////////////////////////////////////////////////////////
//
 
template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options, _Index>::gsSparseMatrix() : Base() { }

template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options, _Index>::gsSparseMatrix(_Index rows, _Index cols) : Base(rows,cols) { }

template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options, _Index>::~gsSparseMatrix() { }

template<typename T, int _Options, typename _Index> inline
void gsSparseMatrix<T, _Options, _Index>::setFrom( gsSparseEntries<T> const & entries) 
        { this->setFromTriplets(entries.begin(),entries.end() ); }

template<typename T, int _Options, typename _Index> inline
gsSparseMatrix<T, _Options,_Index> * gsSparseMatrix<T, _Options, _Index>::clone() const
        { return new gsSparseMatrix(*this); }



}; // namespace gismo
