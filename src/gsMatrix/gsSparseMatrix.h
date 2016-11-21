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

/**
   @brief Class that provides a container for triplets (i,j,value) to be
   filled in a sparse matrix.  

   Constructing a sparse matrix from triplets is much faster than
   inserting directly.  Use gsSparseMatrix().setFrom(gsSparseEntries)
   to pass the triplets to the matrix.

   \tparam T coefficient type
   \ingroup Matrix
*/
template<typename T>
class gsSparseEntries : public std::vector<Eigen::Triplet<T,index_t> >
{
public:
    typedef Eigen::Triplet<T,index_t> Triplet;
    typedef std::vector<Eigen::Triplet<T,index_t> > Base;

    typedef typename Base::iterator       iterator;
public:
    gsSparseEntries() ;

    ~gsSparseEntries();

    inline void add( int i, int j, T value )
    { this->push_back( Triplet(i,j,value) ); }

    inline void addSorted(int i, int j, T value)
    {
        Triplet t(i,j,value);
        iterator pos = std::lower_bound(Base::begin(), Base::end(), t, compTriplet );
        if ( pos == Base::end() || (pos->row() != t.row() && pos->col() !=t.col()) )// If not found
            Base::insert(pos, t);
    }

protected:
    // comparison operator for triplet struct, used to introduce columnwise lex-ordering.
    struct _compTriplet
    {
        bool operator() (const Triplet & left, const Triplet & right)
        {
            return (left.col() < right.col()) ||
                    (left.col() == right.col() && left.row() < right.row()) ;
        }
    } compTriplet ;

};



/** 
    @brief Iterator over the non-zero entries of a sparse matrix

    This class is similar to Eigen::SparseMatrix::InnerIteretor but in
    addition it is default-constructible and assignable.
*/
template<typename T, int _Options, typename _Index>
class gsSparseMatrixIter
{
    typedef Eigen::SparseMatrix<T,_Options,_Index> SparseMatrix;
    static const int IsRowMajor = SparseMatrix::IsRowMajor;

public:
    gsSparseMatrixIter()
    : m_values(NULL), m_indices(NULL), 
      m_end(NULL)   , m_outer(0)
    { }

    gsSparseMatrixIter(const SparseMatrix& mat, const _Index outer)
    : m_outer(outer)
    {
        const _Index oind = mat.outerIndexPtr()[outer];
        m_values  = mat.valuePtr()      + oind;
        m_indices = mat.innerIndexPtr() + oind;
        m_end = mat.isCompressed() ? 
            mat.innerIndexPtr() + mat.outerIndexPtr()[outer+1] : 
            mat.innerIndexPtr() + oind + mat.innerNonZeroPtr()[outer];
    }
    
    inline gsSparseMatrixIter& operator++() 
    { ++m_values; ++m_indices; return *this; }
    
    inline const T& value() const { return *m_values; }
    inline T& valueRef() { return const_cast<T&>(*m_values); }
    
    inline _Index index() const { return *m_indices; }
    inline _Index outer() const { return m_outer; }
    inline _Index row() const { return IsRowMajor ? m_outer : index(); }
    inline _Index col() const { return IsRowMajor ? index() : m_outer; }

    inline operator bool() const { return (m_indices < m_end); }
    
protected:
    const T      * m_values;
    const _Index * m_indices;
    const _Index * m_end;
    _Index m_outer;
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

    typedef gsSparseMatrixIter<T,_Options,_Index> iterator;

    // Type pointing to a block of the sparse matrix
    typedef typename Eigen::Block<Base> Block;

    // Type pointing to a const block of the sparse matrix
    typedef typename Eigen::Block<const Base> constBlock;

    // Type pointing to a block view of the sparse matrix
    typedef gsMatrixBlockView<Base> BlockView;

    /// Shared pointer for gsSparseMatrix
    typedef typename memory::shared<gsSparseMatrix>::ptr Ptr;

    /// Unique pointer for gsSparseMatrix
    typedef typename memory::unique<gsSparseMatrix>::ptr uPtr;

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
    
    /**
       \brief This function returns a smart pointer to the
       matrix. After calling it, the matrix object becomes empty, ie
       the size of the matrix is 0
    */
    Ptr moveToPtr()
    {
        Ptr m(new gsSparseMatrix<T>);
        m->swap(*this);
        return m;
    }
    
    /// \brief Returns an iterator to the first non-zero elemenent of
    /// column \ a outer (or row \a outer if the matrix is RowMajor)
    inline iterator begin(const index_t outer) const { return iterator(*this,outer);}

    void clear()
    {
        this->resize(0,0);
        this->data().squeeze();
    }

    std::pair<index_t,index_t> dim() const 
    { return std::make_pair(this->rows(), this->cols() ); }

    void reservePerColumn(const index_t nz)
    {
        this->setZero();
        this->reserve(gsVector<index_t>::Constant(this->cols(), nz));
    }

    void setFrom( gsSparseEntries<T> const & entries) ;

    inline T   at (_Index i, _Index j ) const { return this->coeff(i,j); }
    inline T & at (_Index i, _Index j ) { return this->coeffRef(i,j); }

    inline T    operator () (_Index i, _Index j ) const { return this->coeff(i,j); }
    inline T  & operator () (_Index i, _Index j ) { return this->coeffRef(i,j); }

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


template<typename T> inline
gsSparseEntries<T>::gsSparseEntries() : Base() { }

template<typename T> inline
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



} // namespace gismo


namespace Eigen { namespace internal {
template<typename T, int _Options, typename _Index>
struct traits<gismo::gsSparseMatrix<T,_Options,_Index> >:
Eigen::internal::traits<Eigen::SparseMatrix<T,_Options,_Index> > { };
} }
